# Variant Calling Algorithm Based On Binomial Distribution

Designed variant calling algorithm in python language. 

## Getting Started

Script expect input file in pileup format described here:

http://samtools.sourceforge.net/pileup.shtml

And outputs only variants in vcf file format on system standard output.

Arguments for running script:

| Argument | Description |
| ------------- |-------------|
|-m, --mpileup 					  | input file in pileup format |
|-t, --threads 					  | number of threads to run in parallel, default is 4 |
|-l, --lines   					  | number of lines to be processed in one chunk, default is 1000 |
|-b, --binomial-distribution-limit | limit in binomial distribution, default is 0.5 |

```
python3 main.py -m test.mpileup                  
```

### Prerequisites

For running script you need python 3.6 version and to install requirements.txt file

```
pip3 install -r requirements.txt
```

## Running the tests

Test can be ran with pytest module

```
python3 -m pytest tests.py
```

## Presentations

Youtube presentation: https://www.youtube.com/watch?v=nufh9DOhOCo&feature=youtu.be


## Authors

* **Nemanja Mikic** 
* **Vidoje Zeljic** 
