# Wrapper script to start processing

# If no input is defined
if [ -z "${INPUT}" ]; then
  export INPUT='/input'
fi

# If no output is defined.
if [ -z "${OUTPUT}" ]; then
 export OUTPUT='/output'
fi

python3 /wrapper_docker.py

#Output log in file system
if test -f "/log.txt"; then
  cp "/log.txt" "$OUTPUT/log.txt"
fi

#Output database in sqlite
if test -f "/processing_db.sqlite"; then
  cp "/db/processing_db.sqlite" "$OUTPUT/processing_db.sqlite"
fi
