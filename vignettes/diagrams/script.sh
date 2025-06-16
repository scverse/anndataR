#!/bin/bash

# Convert mermaid diagrams to different formats
# because RMarkdown doesn't support mermaid diagrams

docker run --rm \
  -w /pwd \
  -v `pwd`:/pwd \
  -u $(id -u):$(id -g) \
  minlag/mermaid-cli \
  -i /pwd/vignettes/diagrams/class_diagram.mmd \
  -o /pwd/vignettes/diagrams/class_diagram.svg
