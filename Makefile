# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# This file specifies the steps to run and their order and allows running them.
# Type `make` for instructions. Type make <command> to execute a command.
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

.DEFAULT_GOAL := help

NAME=$(shell basename `pwd`)
SAMPLES=$(shell ls data)

help:  ## Display help and quit
	@echo Makefile for the $(NAME) package.
	@echo Available commands:
	@grep -E '^[0-9a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m\
		%s\n", $$1, $$2}'

requirements:  ## Install Python requirements
	pip install -r requirements.txt

backup_time:
	echo "Last backup: " `date` >> _backup_time
	chmod 700 _backup_time

_sync:
	rsync --copy-links --progress -r \
	. afr4001@pascal.med.cornell.edu:projects/$(NAME)

sync: _sync backup_time ## [dev] Sync data/code to SCU server

install:
	pip install --use-feature=in-tree-build .

docs:
	cd docs; make html; xdg-open build/html/index.html

clean:
	cd docs; make clean

.PHONY : help \
	requirements \
	sync \
	backup_time \
	_sync \
	sync \
	install \
	docs \
	clean
