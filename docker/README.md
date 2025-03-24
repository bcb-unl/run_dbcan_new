# run_dbcan Docker Image

This Dockerfile builds a container for run_dbcan.
## Build

To build the image:

```bash
docker build -t run_dbcan:v5.0.0 .
```

You can specify a different version using the `BUILD_VERSION` argument:

```bash
docker build --build-arg BUILD_VERSION={your_version} -t run_dbcan:{your_version} .
```

## Usage

Run the container:

```bash
docker run --rm -it run_dbcan:v5.0.0 run_dbcan --help



## Features

- Created based on `micromamba` image
- Support both the `linux/amd64` (Linux and Windows, macOS(Intel)) and `linux/arm64` (macOS(Apple Silicon)).

## License

This docker image is released under the MIT License.
