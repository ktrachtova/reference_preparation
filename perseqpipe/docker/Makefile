include ../config.mk

all: docker_build

docker_build: Dockerfile
	docker build --platform linux/x86_64 -t reference_preparation .

docker_tag:
	docker tag reference_preparation:latest ktrachtok/reference_preparation:latest
	docker tag reference_preparation:latest ktrachtok/reference_preparation:$(REFERENCE_PREPARATION_VERSION)

docker_push:
	docker push --all-tags ktrachtok/reference_preparation

