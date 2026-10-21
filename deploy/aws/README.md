# Running scflow on AWS Batch

> Status: scaffolding. The template and `scflow submit --backend aws-batch` have been
> checked for consistency with `--dry-run`, but have not yet been deployed against a
> real account. Please report problems.

Each job runs on a single EC2 instance, with MPI inside the container. QE/EPW scale
well up to one node for the cell sizes scflow targets, and single-node jobs avoid
the need for EFA networking.

## One-time setup

```bash
# 1. Infrastructure (pick two subnets + a security group from your VPC)
aws cloudformation deploy --stack-name scflow \
  --template-file deploy/aws/batch.yaml --capabilities CAPABILITY_IAM \
  --parameter-overrides VpcSubnets=subnet-aaa,subnet-bbb SecurityGroup=sg-ccc

aws cloudformation describe-stacks --stack-name scflow --query 'Stacks[0].Outputs'

# 2. Push the container image to the ECR repository the stack created
REPO=<RepositoryUri output>
aws ecr get-login-password | docker login --username AWS --password-stdin ${REPO%/*}
docker build -f docker/Dockerfile -t scflow:latest .
docker tag scflow:latest $REPO:latest && docker push $REPO:latest
```

## Submitting

```bash
pip install 'scflow[aws]'
scflow submit compounds/H3S_200GPa.yaml --backend aws-batch \
  --bucket <Bucket> --job-queue <JobQueue> --job-definition <JobDefinition>
```

The job downloads the config from S3 and runs `scflow run` (with the same
pre-flight and validation gate as a local run). It then uploads the inputs and
outputs, the SQLite database and the JSON record to `s3://<Bucket>/scflow/<prefix>/`.
Large binary intermediates (wavefunctions, `.epmatwp`) stay on the instance.

## Cost notes

- Spot instances are on by default. The job definition retries host interruptions,
  but the retry starts from the beginning: the run directory lives on the instance.
  For runs longer than a few hours, set `UseSpot=false`.
- A small cubic cell (Nb, H₃S) finishes in well under an hour on 32 vCPUs. Clathrate
  hydrides on a 4×4×4 q-grid take several hours.
