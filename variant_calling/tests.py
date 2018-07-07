from variant_calling.main import process_line


# test no variant
def test_line_without_variant():
    line = '21	10921797	T	6	....^].^].	?=>?=>'
    result = process_line(0.5, line)
    assert result is None


# if limit is too low and there no other base it will not be variant
def test_line_without_variant_lower_limit():
    line = '21	10921809	C	19	...................	IHHHHIHHHGGH?ECDEAB'
    result = process_line(0, line)
    assert result is None

# test normal case with many difference in reads from ref_base
def test_line_with_variant():
    line = '21	10863087	T	296	G$G$GGggggggggGGggggGGggggGggggggggGGGGgggggGGgggggggGGGggggggggggggggggggGGggggggggGGGggggggggGgggggggggggggggggggggggggggGgggggggggggggggggggggggGGgggggggGggggGggggggggggggGggggggggGggggggggggggggggggggggggGggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg	MJMG@C@C@CCBHjCCFCnGFDFFIE:GCEFEGCIJaGHHGGqkHHIH8DGEnrHIDIHHEIEGI6=<HDFHIeDJIGHH?HrJsHGIDJHHHcICFIIHC5EEJHJIHIHIIIHHIEIJFrICIJIII?JEHIFIIBEC:HJ4HJIEI@JIHJGJJHIdFHJJBIIJIE@Jn;JIHJIIEIJJGICJJIIJIJECJIIJJICJHBEJIIJEEIIIIJJJIFIIF7EIIBJFIII3FIEJIIICIJIJJIIFFCDIGJJ=IIGIIJEJJIJG@JJIIJHDJJIJFJJDJIIIIJA>'
    result = process_line(0.5, line)
    assert result == '21\t10863087\t.\tT\tG\tGT\t1/1'

# if limit is high genotype has more chance to be type 0/1
def test_line_with_variant_lower_limit():
    line = '21	9865961	C	21	.$.TTt,..TT.TT..TT.T,T	DA0<;HHF@@H@?ICA@B@J?'
    result = process_line(0.005, line)
    assert result == '21\t9865961\t.\tC\tT\tGT\t0/1'


# if limit is high genotype has more chance to be type 1/1
def test_line_with_variant_high_limit():
    line = '21	9865961	C	21	.$.TTt,..TT.TT..TT.T,T	DA0<;HHF@@H@?ICA@B@J?'
    result = process_line(0.9, line)
    assert result == '21\t9865961\t.\tC\tT\tGT\t1/1'


# if limit is too low every like here if there was only one base is different
# from ref_base it will be variant
def test_binomial_distribution_line():
    line = '21	9865789	G	36	.$.,...,,.a...,..,,....,......,,,....	B>ElFFEHFABGGEdFHBE@FFCFGFG05CE=EFF@'
    result = process_line(0, line)
    assert result == '21\t9865789\t.\tG\tA\tGT\t0/1'
