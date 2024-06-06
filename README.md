# changeseq-snakemake
 This repository is a Snakemake-based workflow designed to streamline gene editing detection processes, reimplementing the original Changseq repository for enhanced efficiency and manageability.

## Installation

Ensure you have Python 3 installed and the following dependencies:

- Snakemake
- required_rules_module

You can install the dependencies for this project using the following command:

```bash
pip install -r requirements.txt
```

Usage Instructions
Clone this repository:
```
git clone https://github.com/your_username/changeseq-snakemake.git
cd changeseq-snakemake
```
Prepare your input data:

Ensure you have sufficient input data, such as sequencing files and a reference genome.

Run the analysis:

Start the analysis by running the following command in the root directory of the project:
```
snakemake --use-conda
```
For more customized options, please refer to the official Snakemake documentation.

Contribution Guidelines
We welcome all forms of contribution, including bug reports, suggestions, and improvement proposals. Please use the GitHub issue tracker or submit a PR to make a contribution.

License
This repository is licensed under the MIT License. Please refer to the LICENSE file for detailed information.
