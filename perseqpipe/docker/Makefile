include ../config.mk

IMAGE_NAME=reference_preparation
DOCKER_USER=ktrachtok

docker_build_x86_64: Dockerfile
	docker buildx build \
		--platform linux/x86_64 \
		-t $(IMAGE_NAME):x86_64-$(REFERENCE_PREPARATION_VERSION) \
		-f Dockerfile ..

# Tag image
docker_tag_x86_64:
	docker tag $(IMAGE_NAME):x86_64-$(REFERENCE_PREPARATION_VERSION) $(DOCKER_USER)/$(IMAGE_NAME):x86_64-$(REFERENCE_PREPARATION_VERSION)

# Push image
docker_push_x86_64:
	docker push $(DOCKER_USER)/$(IMAGE_NAME):x86_64-$(REFERENCE_PREPARATION_VERSION)

# Optional: Build all and push all
docker_all: docker_build_x86_64 docker_tag_x86_64 docker_push_x86_64