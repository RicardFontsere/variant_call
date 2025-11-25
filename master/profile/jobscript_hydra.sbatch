#!/bin/bash

# VUB-HPC specific customizations - do not modify!

OS_SUBMIT=$VSC_OS_LOCAL
ARCH_SUBMIT=$VSC_ARCH_LOCAL
SUFFIX_SUBMIT=$VSC_ARCH_SUFFIX


# workaround to make --export=ALL work in Hydra on different partitions
unset VSCPROFILELOADED VSC_OS_LOCAL VSC_ARCH_LOCAL VSC_ARCH_SUFFIX
source /etc/profile.d/vsc.sh
module update

# workaround to make sure we use the correct python executable in compute nodes
envpython=${{EBROOTPYTHON/$OS_SUBMIT\/$ARCH_SUBMIT$SUFFIX_SUBMIT/$VSC_OS_LOCAL\/$VSC_ARCH_LOCAL$VSC_ARCH_SUFFIX}}/bin/python
eval "$envpython() {{ python \"\$@\"; }}"


#########################################################################################

# properties = {properties}
{exec_job}


