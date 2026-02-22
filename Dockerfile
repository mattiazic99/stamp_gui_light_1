# ==========================
# 1) Immagine base Python
# ==========================
FROM python:3.11-slim

# ==========================
# 2) Directory di lavoro dentro il container
# Tutto il tuo progetto sarà copiato in /app
# ==========================
WORKDIR /app

# ==========================
# 3) Copia e installa i requirements
# (più veloce e più pulito)
# ==========================
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# ==========================
# 4) Copia tutto il progetto
# ==========================
COPY . .

# ==========================
# 5) Variabili Streamlit
# ==========================
ENV STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
ENV STREAMLIT_SERVER_ENABLECORS=true
ENV STREAMLIT_SERVER_ENABLEXSRS=false

# ==========================
# 6) Espone la porta di Streamlit
# ==========================
EXPOSE 8501

# ==========================
# 7) Comando di avvio
# ==========================
CMD ["streamlit", "run", "main.py", "--server.port=8501", "--server.address=0.0.0.0"]
