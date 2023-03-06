# DnaStorage

## Usage

Create a venv (recommended)

```console
git clone https://github.com/InbalPreuss/DnaStorage.git
pip install -r requirements.txt
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


