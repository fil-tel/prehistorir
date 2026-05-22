##################################################################
# Docker shortcuts
#

# check the OS type to assign appropriate Docker image tag
ifeq ($(shell uname -s), Darwin)
    PLATFORM ?= arm64
else
    PLATFORM ?= amd64
endif

IMAGE := $(shell basename $(shell pwd)):$(PLATFORM)
CONTAINER := $(shell basename $(shell pwd))

# if present, extract GitHub access token
TOKEN := $(shell if [[ -f ~/.GITHUB_PAT ]]; then more ~/.GITHUB_PAT; else echo ""; fi)

.PHONY: rstudio bash r docker-build docker-push docker-pull local-webapp remote-webapp

bash:
	docker run --rm -ti -v $(shell pwd):/project -w /project --name $(CONTAINER) $(IMAGE) bash

r:
	docker run --rm -ti -v $(shell pwd):/project -w /project --name $(CONTAINER) $(IMAGE) R

attach:
	container_id=`docker ps | awk -v name="$$(basename "$$PWD")" '$$2 ~ name {print $$1}'`; \
	docker exec -it $$container_id /bin/bash

rstudio:
ifndef PORT
	$(error PORT variable must be set explicitly)
endif
	docker run --rm -ti -p $(PORT):8787 -e RUNROOTLESS=true -e DISABLE_AUTH=true -v $(shell pwd):/project -w /project --name $(CONTAINER) $(IMAGE)

docker-build:
	docker build --build-arg GITHUB_PAT=$(TOKEN) -t $(IMAGE) .

docker-build-clean:
	@echo $(IMAGE)
	docker build --no-cache --build-arg GITHUB_PAT=$(TOKEN) -t $(IMAGE) .

docker-push:
	docker push $(IMAGE)

docker-pull:
	docker pull $(IMAGE)

docker-stop:
	docker stop $(CONTAINER) || true

local-webapp:
ifndef PORT
	$(error PORT variable must be set explicitly)
endif
	firefox http://localhost:$(PORT) &

remote-webapp:
ifndef SERVER
	$(error SERVER variable must be set to start a web app)
endif
ifndef PORT
	$(error PORT variable must be set explicitly)
endif
	@if [[ "$(SERVER)" != "localhost" ]]; then \
	    PID=$$(lsof -ti:$(PORT)); \
	    if [[ -n "$$PID" ]]; then \
		kill -9 $$PID; \
	    fi; \
	    autossh -M 0 -f -N -L localhost:$(PORT):localhost:$(PORT) $(SERVER) || { \
		echo "SSH connection failed. Exiting."; \
		exit 1; \
	    }; \
	fi; \
	firefox http://localhost:$(PORT) &

port-forward:
ifndef SERVER
	$(error SERVER variable must be set to start a web app)
endif
ifndef PORT
	$(error PORT variable must be set explicitly)
endif
	@if [[ "$(SERVER)" != "localhost" ]]; then \
	    PID=$$(lsof -ti:$(PORT)); \
	    if [[ -n "$$PID" ]]; then \
		kill -9 $$PID; \
	    fi; \
	    autossh -M 0 -f -N -L localhost:$(PORT):localhost:$(PORT) $(SERVER) || { \
		echo "SSH connection failed. Exiting."; \
		exit 1; \
	    }; \
	fi
