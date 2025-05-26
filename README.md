# scRNA_seq_Pipeline

# Οδηγίες Εκτέλεσης της Εφαρμογής (μόνο για Windows)

Ακολουθήστε τα παρακάτω βήματα για να τρέξετε την εφαρμογή μέσω Docker:

---

## 🐳 1. Άνοιγμα Docker Desktop

> Ανοίξτε την εφαρμογή **Docker Desktop** και βεβαιωθείτε ότι τρέχει σωστά.

---

## 💻 2. Άνοιγμα Terminal

> Μεταβείτε στον φάκελο `bioml_streamlit_app` μέσω **PowerShell** ή **Command Prompt**.

---

## ⚙️ 3. Δημιουργία Docker Image

```bash
docker build -t scrna-app .
```

Αυτό δημιουργεί ένα Docker image με όνομα `scrna-app`.

---

## 🚀 4. Εκτέλεση της Εφαρμογής

```bash
docker run -p 8501:8501 scrna-app
```

Αυτό θα ξεκινήσει τον Streamlit server μέσα στο container.

---

## 🌐 5. Άνοιγμα στο Browser

> Ανοίξτε τον browser σας και επισκεφθείτε τη διεύθυνση:

```
http://localhost:8501/
```

---

Εφόσον όλα έχουν ρυθμιστεί σωστά, η εφαρμογή θα είναι πλήρως λειτουργική μέσω του Docker.
