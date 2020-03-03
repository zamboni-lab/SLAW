###We mount sauer1 if any or the input or output is running on sauer 1
###If no input is defined
if [ -z "${INPUT}" ]; then
  export INPUT='/input'
fi

###If no output is defined.
if [ -z "${OUTPUT}" ]; then
 export OUTPUT='/output'
fi

##Linux command to get available ram evnetuall


#If the output is one sauer1 we try to mount the disk
if [[ "$OUTPUT" == /sauer1*  ]] || [[ "$INPUT" == /sauer1*  ]] ;
then
  mkdir sauer1
  ##In this case we mount sauer 1
  if [[ -z "${PASSWORD}" ]]; then
    #If the password is defined we mount it
    # echo "mount //nas22.ethz.ch/biol_imsb_sauer_1 sauer1 -nobrl -rw -o domain=d.ethz.ch,username=$USERNAME,vers=3.0"
    mount //nas22.ethz.ch/biol_imsb_sauer_1 sauer1 -nobrl -rw -o domain=d.ethz.ch,username=$USERNAME,vers=3.0
  else
    # echo "mount //nas22.ethz.ch/biol_imsb_sauer_1 sauer1 -nobrl -rw -o domain=d.ethz.ch,username=$USERNAME,vers=3.0,password=$PASSWORD"
    mount //nas22.ethz.ch/biol_imsb_sauer_1 sauer1 -nobrl -rw -o domain=d.ethz.ch,username=$USERNAME,vers=3.0,password=$PASSWORD
  fi

  ## WE ehck that suaer1 has been mounted ocrrectly
  if mountpoint -q -- "sauer1"; then
	echo "sauer1 mounted correctly"
  else
	echo "sauer1 does not seems to have been mounted correctly, please check your password and username ?"
  fi

fi

if [ "$(ls -A $OUTPUT)" ]; then
     echo "Destination directory is not empty."
fi

python3 wrapper_docker.py
if test -f "/log.txt"; then
  cp "/log.txt" "$OUTPUT/log.txt"
fi

if test -f "/processing_db.sqlite"; then
  cp "/db/processing_db.sqlite" "$OUTPUT/processing_db.sqlite"
fi
