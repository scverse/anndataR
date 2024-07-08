#!/bin/bash

# Convert mermaid diagrams to different formats
# because RMarkdown doesn't support mermaid diagrams

docker run --rm -w /pwd -u `id -u`:`id -g` \
  -v `pwd`:/pwd minlag/mermaid-cli \
  -i vignettes/diagrams/class_diagram.mmd -o vignettes/diagrams/class_diagram.svg
