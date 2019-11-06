# Installation

```bash
git clone https://github.com/moozeq/DD_SCAFFOLDS.git

cd DD_SCAFFOLDS
pip3 install -r requirements.txt
```

# Usage
```bash
./scaffolds.py <enzyme target> <inhibitors comparing method>
```

## Example
Make AMCase scaffolds with RINGS_WITH_LINKERS_1
```bash
./scaffolds.py amcase RINGS_WITH_LINKERS_1
```

Scaffolds should be displayed on standard output, in folder there should be file **amcase** with all inhibitors in format "SMILE ID"