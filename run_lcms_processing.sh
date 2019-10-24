###We mount sauer1 if any or the input or output is running on sauer 1
###If no input is defined
if [ -z "${INPUT}" ]; then
  export INPUT='/input'
fi

###If no output is defined.
if [ -z "${OUTPUT}" ]; then
 export OUTPUT='/output'
fi

#If the output is one sauer1 we try to mount the disk
if [[ "$OUTPUT" == /sauer1*  ]] || [[ "$INPUT" == /sauer1*  ]] ;
then
  mkdir sauer1
  ##In this case we mount sauer 1
  if [[ -z "${PASSWORD}" ]]; then
    #If the password is defined we mount it
    echo "mount //nas22.ethz.ch/biol_imsb_sauer_1 sauer1 -nobrl -rw -o domain=d.ethz.ch,username=$USERNAME,vers=3.0"
    mount //nas22.ethz.ch/biol_imsb_sauer_1 sauer1 -nobrl -rw -o domain=d.ethz.ch,username=$USERNAME,vers=3.0
  else
    echo "mount //nas22.ethz.ch/biol_imsb_sauer_1 sauer1 -nobrl -rw -o domain=d.ethz.ch,username=$USERNAME,vers=3.0,password=$PASSWORD"
    mount //nas22.ethz.ch/biol_imsb_sauer_1 sauer1 -nobrl -rw -o domain=d.ethz.ch,username=$USERNAME,vers=3.0,password=$PASSWORD
  fi
fi

if [ "$(ls -A $OUTPUT)" ]; then
     echo "Destination directory is not empty."
fi

python3 wrapper_docker.py
