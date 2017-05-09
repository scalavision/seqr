#!/usr/bin/env bash

DATABASE_NAME=$1
BACKUP_FILE=$2

if [ -z $DATABASE_NAME ]; then
    echo Database name not specified
    exit -1
fi

if [ -z $BACKUP_FILE ]; then
    echo Backup file name not specified
    exit -1
fi

psql -U postgres $DATABASE_NAME < <(gunzip -c $BACKUP_FILE)  # load a .txt.gz file generated by pg_dump
