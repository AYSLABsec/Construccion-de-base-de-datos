# Indexador y renombrador de secuencias de Bacillus

Este repositorio facilita el procesamiento de un archivo FASTA con todas las secuencias de Bacillus, permitiendo extraer subconjuntos de secuencias, generar un índice en formato CSV y renombrar los encabezados de las secuencias de forma automatizada.

---

## Estructura del repositorio

- **getSeqByListID.py**\
  Script que recibe un archivo FASTA y un listado de IDs (uno por línea) y genera un FASTA con las secuencias correspondientes. fileciteturn0file0

- **extraxt-seq.py**\
  Función alternativa para extraer secuencias de un FASTA a partir de un archivo de texto con IDs, escribiendo en `extracted_sequences.fasta`. fileciteturn0file1

- **rename.py**\
  Renombra los encabezados de un FASTA usando un archivo CSV de mapeo (ID	Descripción) e inserta el nombre de gen indicado por el usuario. fileciteturn0file2

---

## Requisitos previos

- Python 3.x
- Biopython (`pip install biopython`)

---

## Instalación

```bash
# Clonar el repositorio
git clone <URL_DEL_REPOSITORIO>
cd <NOMBRE_REPOSITORIO>

# Instalar dependencias
pip install -r requirements.txt  # Si no existe, instala Biopython manualmente:
pip install biopython
```

---

## Flujo de trabajo

### 1. Preparar el archivo FASTA con las secuencias completas

Reúne todas las secuencias de Bacillus en un único archivo, por ejemplo `all_bacillus.fasta`.

### 2. Generar un índice CSV

Crea un archivo `index.csv` con dos columnas separadas por tabulador (o coma):

```
sequence_id	descripcion
Bac123	"Proteína X de Bacillus subtilis"
Bac456	"Enzima Y de Bacillus anthracis"
...
```

> **Tip:** puedes extraer rápidamente los IDs originales con:
>
> ```bash
> ```

grep '^>' all\_bacillus.fasta | sed 's/>//' | cut -f1 -d '|' > ids.txt

# Luego editar manualmente o con Excel para completar descripciones.

````
> Asegúrate de guardar el archivo como CSV de texto.

### 3. Extraer subconjuntos de secuencias (opcional)

- Con **getSeqByListID.py**:
  ```bash
  python getSeqByListID.py all_bacillus.fasta ids.txt
  # Salida: all_bacillus_genes.fasta
````

- Con **extraxt-seq.py**:
  ```bash
  python extraxt-seq.py all_bacillus.fasta ids.txt
  # Salida: extracted_sequences.fasta
  ```

### 4. Renombrar encabezados

```bash
python rename.py -i extracted_sequences.fasta -m index.csv -o renamed.fasta
```

El script pedirá interactivamente el **nombre del gen** a incluir en los encabezados. El archivo resultante `renamed.fasta` tendrá encabezados con el formato:

```
>sequence_id\tdescripcion\tgen_name
SECUENCIA...
```

---

## Ejemplo de uso

1. Extraer IDs de interés:
   ```bash
   grep '^>' all_bacillus.fasta | sed 's/>//' | cut -f1 -d '|' > ids_of_interest.txt
   ```
2. Ejecutar extracción:
   ```bash
   python getSeqByListID.py all_bacillus.fasta ids_of_interest.txt
   ```
3. Crear o ajustar `index.csv` con las descripciones.
4. Renombrar secuencias:
   ```bash
   python rename.py -i all_bacillus_genes.fasta -m index.csv -o final_sequences.fasta
   ```

---

## Contribuciones

¡Las contribuciones son bienvenidas! Para reportar errores o sugerir mejoras, abre un *issue* o envía un *pull request*.

---

## Licencia

Este proyecto está licenciado bajo la licencia MIT. Consulta el archivo [LICENSE](LICENSE) para más detalles.

