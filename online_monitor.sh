#!/usr/bin/env bash
## Author V. Balagura, balagura@cern.ch (19.11.2012)


export ONLINE_MONITOR_DIR=`dirname "$0"`
export R_PROFILE_USER=${ONLINE_MONITOR_DIR}/online_monitor.R

TEMP=`getopt -o hp: --long help,pedestal-suppression: -n $0 -- "$@"`

if [ $? != 0 ] ; then echo "Usage: $0 [-h|--help] [(-p|--pedestal-suppression) 0.10] [file.raw]" ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

while true ; do
        case "$1" in
                -h|--help) echo "Usage: $0 [-h|--help] [(-p|--pedestal-suppression) 0.10] [file.raw]" ; shift ;;
                -p|--pedestal-suppression) echo "Not triggered hits are prescaled by $2"; export ONLINE_MONITOR_PEDESTAL_SUPPRESSION=$2 ; shift 2 ;;
                --) shift ; break ;;
                *) echo "Internal error!" ; exit 1 ;;
        esac
done
export ONLINE_MONITOR_FILE=$1

R --no-save --no-restore --quiet

# Another way to do the same (but without need for
# for (pack in options("defaultPackages")[[1]]) require(pack, character.only = TRUE)
# at the top of online_monitor.R):
#
# (
# cat <<'EOF'
# .First <- function() {
#   for (pack in options("defaultPackages")[[1]]) require(pack, character.only = TRUE) # load default libs
#   source("online_monitor.R", echo=TRUE)
# }
# EOF
# ) > /home/balagura/c/ecal/online_monitor/.Rprofile
# R --no-save
# rm -f /home/balagura/c/ecal/online_monitor/.Rprofile
