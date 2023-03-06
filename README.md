# DnaStorage

## Usage

Create a venv (recommended)

```console
git clone https://github.com/InbalPreuss/DnaStorage.git
pip install -r requirements.txt
```

inside DnaStorage directory, next to dna_storage directory create data/testing/input_text.dna with some text inside.

All the outputs will be written next to it.

then just run:
```console
python main.py
```

## Test
```console
python -m pytest
```
Create graphs:
In test_dna.py run test_number_of_oligos_per_barcode()

## Profiling
in test_dna.py make sure that the function in
```console
if __name__ == '__main__':
```
is:
```console
code_profiling()
```
Then, run:

linux:
```console
python -m cProfile -s time test_dna.py > temp_file && head -n 10 temp_file
```
Windows: 
```console
python -m cProfile -s time test_dna.py > profiling_data_no_synthsis_1_KB.txt
```

##Run simulation experiment:
1. Create simulation data. The data will be in data/testing.
```console
nohup python3 -m tests.distributed &
```
2. Run plots on the data created. The data will be in .
When finished run:
```console
nohup python3 -m dna_storage.plots &
```


