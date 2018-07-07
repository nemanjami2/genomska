import io
import argparse
import functools
from collections import defaultdict
from multiprocessing import Pool
from scipy.stats import binom

bases = {'A', 'C', 'T', 'G'}

same_as_ref = {'.', ','}

allowed_bases = bases.union(same_as_ref)

special_characters = {'+', '-'}

error_probabilities = {}


class ParsedBase:

    def __init__(self, count=0, quality=0):
        self.count = count
        self.quality = quality


def process_line(binomial_limit, line):
    """
    Calls SNP for single line
    """
    data = line.split()
    parsed_data = defaultdict(ParsedBase)

    # data[0]   [1] [2]     [3]         [4]         [5]
    # reference pos refBase coverage    read bases  baseQualities

    ref_base = data[2]
    qualities = data[5]
    index = 0
    special_character = None
    was_special_case = False

    # Convert bases to upper
    pileup_bases = data[4].upper()
    old_base = None
    for read_base in pileup_bases:

        # special case ignore
        if was_special_case:
            was_special_case = False
            continue

        # ^ marks quality of read
        if read_base == '^':
            was_special_case = True
            continue

        # if it's fist base
        if not old_base:
            if read_base in same_as_ref:
                read_base = ref_base
            old_base = read_base
            continue

        # continue to parse insertion or deletion
        if special_character:
            special_character = special_character + read_base
            num_to_calc = int(special_character[2])
            num_index = 3
            while (
                len(special_character) > num_index
                and special_character[num_index].isdigit()
            ):
                num_to_calc = num_to_calc * 10 + int(
                    special_character[num_index]
                )
                num_index = num_index + 1

            if len(special_character) == (
                2 + len(str(num_to_calc)) + num_to_calc
            ):
                old_base = special_character
                special_character = None
            continue

        # special character is + and - it means insertions or deletion in that
        # position
        if read_base in special_characters:
            special_character = old_base + read_base
            continue

        # if read_base is . or , it is ref_base
        if read_base in same_as_ref:
            read_base = ref_base

        # skip bases with p(error)=1
        if qualities[index] == '!':
            index += 1
            continue

        # skip letters that are not bases, e.g. start or end of read markers
        if read_base not in allowed_bases:
            continue

        # save info for read base in previous iteration
        old_data = parsed_data[
            old_base
        ]
        old_data.count += 1
        old_data.quality += ord(qualities[index])
        index += 1
        parsed_data[old_base] = old_data
        old_base = read_base

    # save for last read base
    old_data = parsed_data[
        old_base
    ]
    old_data.count += 1
    old_data.quality += ord(qualities[index])
    index += 1
    parsed_data[old_base] = old_data

    # sort all bases by number of repetition
    sorted_by_value = sorted(
        parsed_data.items(), key=lambda kv: kv[1].count, reverse=True
    )

    # base with the most repetition
    call = sorted_by_value[0][0]

    kind_of_variant = '1/1'

    # if there is more than one baes
    if len(sorted_by_value) > 1:
        kind_of_variant = '0/1'
        # calculate count for binomial distribution
        n = sorted_by_value[0][1].count + sorted_by_value[1][1].count
        value_to_check = sorted_by_value[0]
        second_value = sorted_by_value[1]
        p1 = sorted_by_value[1][1].quality/sorted_by_value[1][1].count - 33
        p0 = sorted_by_value[0][1].quality/sorted_by_value[0][1].count - 33

        # if the most repeatable is ref_base swap values
        if value_to_check[0] == ref_base:
            value_to_check, second_value = second_value, value_to_check
            p1, p0 = p0, p1


        # calculate k and p for binomial distribution
        k = value_to_check[1].count
        p = p0/p1

        # calculate limit for binomial distribution
        limit = binom.pmf(int(n*p), n, p) * binomial_limit
        result_call = binom.pmf(k, n, p)

        # check if there was variant
        if result_call >= limit or k >= n/2:
            call = value_to_check[0]
            if result_call < limit:
                if second_value[0] != ref_base:
                    kind_of_variant = '1/2'
                else:
                    kind_of_variant = '1/1'
        else:
            call = second_value[0]
            if call != ref_base:
                kind_of_variant = '1/2'

    if call != ref_base and call != '*':
        if call.find('+') > 0:
            call = call.replace('+', '')
            call = ''.join([i for i in call if not i.isdigit()])
        elif call.find('-') > 0:
            call = call.replace('-', '')
            call = ''.join([i for i in call if not i.isdigit()])
            if ref_base == call[0]:
                ref_base, call = call, ref_base
        return '\t'.join(
            [data[0], data[1], '.', ref_base, call, 'GT',
             kind_of_variant]
        )


class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        return self.start <= other <= self.end


class Config:

    def __init__(self):
        self.input_file = ''
        self.threads = 4
        self.lines = 1000
        self.binomial_limit = 0.5

    def set_data_from_dict(self, arguments: dict):
        self.lines = arguments.get('lines')
        self.threads = arguments.get('threads')
        self.input_file = arguments.get('mplieup')
        self.binomial_limit = arguments.get(
            'binomial_distribution_limit'
        )
        return self


def command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-m', '--mplieup', required=True,
        type=str, help='Input file in pileup format'
    )
    parser.add_argument(
        "-t", "--threads", default=4,
        type=int, help='Number of threads'
    )
    parser.add_argument(
        "-l", "--lines", default=1000,
        type=int, help='How many lines in bucket for read'
    )
    parser.add_argument(
        "-b", "--binomial-distribution-limit", default=0.5,
        type=float, help='Limit for binomial distribution,',
        choices=[Range(0.0, 1.0)]
    )
    args = vars(parser.parse_args())
    return Config().set_data_from_dict(args)


def print_header():
    print('##fileformat=VCFv4.2')
    print('#CHROM\tPOS\tID\tREF\tALT\tFORMAT\tHCC1143BL')


if __name__ == "__main__":
    config = command_line()
    pool = Pool(processes=config.threads)
    process_line_with_config = functools.partial(
        process_line, config.binomial_limit
    )
    try:
        with io.open(config.input_file, 'r') as file:
            print_header()
            queue = []
            for line in file.readlines():
                if len(queue) < config.lines:
                    queue.append(line)
                else:
                    result = pool.map(process_line_with_config, queue)
                    for item in result:
                        if item:
                            print(item)
                    queue = []

            for item in pool.map(process_line_with_config, queue):
                if item:
                    print(item)
    except FileNotFoundError:
        print(
            'File with file path "{file_name}" doesn\'t exists'.format(
                file_name=config.input_file
            )
        )
