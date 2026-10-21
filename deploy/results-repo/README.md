# scflow results

Superconducting Tc calculations submitted by people running
[scflow](https://github.com/Vdamarla7/Superconductivity).

Each file in `records/<formula>/<pressure>GPa/` is one run. It contains:
- the structure and all calculation settings;
- λ, ω_log and the Tc estimates;
- the outcome of every scflow check;
- code versions and hashes of the input and output files.

It contains no hostnames or file paths.

**These results are self-reported.** The tier (`invalid` / `preliminary` / `converged` /
`benchmark`) was computed by scflow on the submitter's machine. A record is only marked verified
after someone reproduces the run and its id is added to `verified.txt`.

Browse the table at the GitHub Pages site for this repo.

## Submitting

scflow submits finished runs automatically when you're logged in to GitHub (`gh auth login` or
`GITHUB_TOKEN`). It opens a pull request from your fork. Queued results can be sent with
`scflow share`. CI checks each record against `schema/record.v1.json` before it is merged.

## License

The data in `records/` is licensed under CC BY 4.0
(https://creativecommons.org/licenses/by/4.0/). The scripts are MIT.
