#!/bin/bash

if [[ "$TESTENV" == "build" ]]
then
	echo "$DOCKER_PASSWORD" | docker login -u "$DOCKER_USERNAME" --password-stdin
	docker push tonyyzy/gc_analysis
fi
