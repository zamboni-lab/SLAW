docker build -t lcms_workflow_zamboni .

docker image rm 4f7c30049e8b

docker run -v /var/run/docker.sock:/var/run/docker.sock -v E:/output_docker:/output --privileged -t --rm quay.io/singularity/docker2singularity:v2.6 lcms_workflow_zamboni
