FROM python:3.12-slim

WORKDIR /app

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

# SECURITY: AUDIT_SECRET_KEY must be provided at runtime, not baked into image
# Generate a secure key: python -c "import secrets; print(secrets.token_hex(32))"
# Example: docker run -e AUDIT_SECRET_KEY=<your-secure-key> her2-ihc-fish-algorithm

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt* ./
RUN pip install --no-cache-dir fastapi uvicorn pydantic pytest

COPY . .

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

CMD ["python", "cli.py", "serve", "--host", "0.0.0.0", "--port", "8000"]
