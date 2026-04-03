# seqair development container with Claude Code
#
# Usage:
#   docker build -t seqair-dev .
#
#   # tmpfs target/ — fastest builds, cold start on every run (RAM-backed, lost on exit)
#   docker run -it --rm \
#     -v "$HOME/.claude:/home/dev/.claude" \
#     -v "$(pwd):/workspace" \
#     --tmpfs /workspace/target:rw,size=8g \
#     -w /workspace \
#     seqair-dev
#
#   # Named volume target/ — fast builds, warm cache persists across runs
#   docker run -it --rm \
#     -v "$HOME/.claude:/home/dev/.claude" \
#     -v "$(pwd):/workspace" \
#     -v seqair-target:/workspace/target \
#     -w /workspace \
#     seqair-dev
#
# Then inside: claude --dangerously-skip-permissions

FROM debian:bookworm-slim

# System packages
# - build-essential, cmake, autoconf, automake, pkg-config: C build toolchain
#   (rust-htslib static feature compiles htslib from C source)
# - zlib1g-dev, libbz2-dev, liblzma-dev: compression libs linked by htslib + bzip2/xz2 crates
# - libcurl4-openssl-dev, libssl-dev: htslib remote I/O (needed even with static feature)
# - git: cargo fetches the rust-htslib git dependency at build time
# - curl, ca-certificates: download rustup and Node setup script
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    autoconf \
    automake \
    libtool \
    pkg-config \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    git \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Node.js 22 (LTS) — required by Claude Code
RUN curl -fsSL https://deb.nodesource.com/setup_22.x | bash - \
    && apt-get install -y --no-install-recommends nodejs \
    && rm -rf /var/lib/apt/lists/*

# Claude Code
RUN npm install -g @anthropic-ai/claude-code

# Rust installed system-wide so the dev user can use it without rustup in PATH tricks.
# CARGO_HOME/RUSTUP_HOME under /usr/local — world-readable, cargo bin on global PATH.
ENV RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:$PATH
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
    | sh -s -- -y --no-modify-path --default-toolchain stable --profile minimal \
    && rustup component add clippy rustfmt \
    # Allow the non-root user to write into cargo's registry/cache
    && chmod -R a+w /usr/local/cargo

# Non-root user — required for Claude Code --dangerously-skip-permissions (yolo mode)
RUN useradd -m -s /bin/bash dev
USER dev

# Pre-warm the Cargo registry + compile workspace dependencies so that
# incremental rebuilds when mounting source are fast.
# We copy only the manifest files first (Docker layer cache trick).
WORKDIR /home/dev/dep-cache

COPY --chown=dev Cargo.toml Cargo.lock ./
COPY --chown=dev crates/seqair/Cargo.toml crates/seqair/Cargo.toml
COPY --chown=dev crates/seqair-types/Cargo.toml crates/seqair-types/Cargo.toml

# Create stub lib files so cargo can resolve and compile deps
RUN mkdir -p crates/seqair/src crates/seqair-types/src \
    && echo 'fn main() {}' > crates/seqair/src/main.rs \
    && touch crates/seqair/src/lib.rs crates/seqair-types/src/lib.rs \
    && cargo fetch \
    && cargo build --workspace 2>/dev/null || true \
    && cargo test --workspace --no-run 2>/dev/null || true

WORKDIR /workspace

CMD ["bash"]
