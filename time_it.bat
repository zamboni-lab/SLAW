$sw = [Diagnostics.Stopwatch]::StartNew()
docker run -it --cap-add=SYS_ADMIN --cap-add=DAC_READ_SEARCH --privileged -e POLARITY=positive -e INPUT=/sauer1/users/Alexis/examples_lcms_workflow/input -e PASSWORD=wyx5z5r9Milena90 -e OUTPUT=/sauer1/users/Alexis/examples_lcms_workflow/output -e USERNAME=dalexis adelabriere/lcms_workflow_zamboni:stable
$sw.Stop()
$sw.Elapsed


$sw1 = [Diagnostics.Stopwatch]::StartNew()
docker run -it --cap-add=SYS_ADMIN --cap-add=DAC_READ_SEARCH --privileged -e POLARITY=positive -e INPUT=/sauer1/users/Alexis/examples_lcms_workflow/input -e PASSWORD=wyx5z5r9Milena90 -e OUTPUT=/sauer1/users/Alexis/examples_lcms_workflow/output -e USERNAME=dalexis lcms_workflow_zamboni
$sw1.Stop()
$sw1.Elapsed
