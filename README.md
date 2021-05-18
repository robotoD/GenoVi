# GenomeVisualizer

Genome visualizer software project

## Execution:

### create_kar.py

```python
python create_kar.py -i <input genbank file path> -o <output kar file path>
```

Example with testing data: 
```python
python create_kar.py -i input_test/PROKKA_01232018.gbk -o test.kar
``` 
## Meetings journal

### Reunión 06/01

* GBK a fna, GBK a faa, GBK a gff.
* Mirar proyecto [SnakemakeCircos](https://bitbucket.org/mroachawri/snakemakecircos/src/master/). Ver cómo se podría adaptar.
* Ver cómo mejorar la sintaxis.
* Ver qué cosas se puede agregar (Cog).
* Hacer grupo de Telegram

### Reunión 04/05

* Generalmente los archivos están en GBK y usar este archivo es lo mejor ya que te da toda la información que otros archivos te dan por partes (fna, faa y gff).
* Nos importa graficar contenido GC, simetría GC, secuencias codificantes que encontró (verde una hebra y azul la otra) y tamaño del genoma (cada fragmento representa un contig).
* Opcionalmente si el usuario quiere hacer la clasificación COG (que a veces viene en GBK), se necesitan todas las secuencias de nucleótidos, por lo que ahí sale más fácil que el usuario ingrese el archivo fna. Para obtener la anotación se puede correr [DFAST](https://dfast.ddbj.nig.ac.jp/) con el archivo fna. Esto produce un nuevo GBK que sí tiene información sobre COG categories.
* (&) chr.kar es un archivo input para CIRCOS. Este archivo se crea tomando el valor indicado en el campo SOURCE del archivo gbk. La primera línea contiene el tamaño total del genoma. "band" simboliza cada uno de los contigs. Un ejemplo de SOURCE para dos contigs consecutivos: SOURCE_1 "1...N" y SOURCE_2 "1...M". Para ponerlo en chr.kar se debe tomar "0...N" y "N+1...N+M".
* Se podría también graficar los repeat_regions, tRNA, rRNA. Además de destacar ciertas regiones que el usuario puede elegir.

#### Tareas
* Andrea: Tener (&) listo para la próxima reunión.
* Roberto: Revisar lo que hace DFAST. Ver cómo se anotan los COGs. Preguntar a Vicente si se quiere unir al proyecto. (https://github.com/transcript/COG; Script de python solo requiere descargar la base de datos del NCBI de COGs, crear la base de datos de COG con DIAMOND, y luego buscar con DIAMOND en el .faa. Entrega una tabla con 3 columnas, la tercera es el código COG). Vicente esta abierto a participar.
* Andrés: Averiguar qué información se debe extraer del GBK para hacer los archivos input para colorear secciones custom.


### Reunión 18/05

* Prueba de primera versión de create_kar.py
* De create_kar modificar: El archivo .kar final debe quedar en el mismo orden de como aparecen los contigs en el .gbk. El orden de las bandas se cuenta desde 0 de nuevo. Agregar prefijos de argumentos, verificadores y comentarios a create_kar.
* Siguiente archivo a crear es el que define las flechas de las secuencias codificantes (CDSs_positive.txt). En el archivo gbk, en cada contig aparece el campo CDS, el location indica donde inicia y termina la secuencia. Si aparece complement(2..444) significa que la secuencia es negativa. Si aparece 2..444 la secuencia es positiva.
* A la locación de todos los CDS de un contig se les debe sumar el tamaño total del contig anterior --El orden da igual--. Por ejemplo, si dice el contig anterior dice "source 1..100" y el primer CDS del contig siguiente tiene location "4..10", debe quedar como "104..110".
* La última columna del CDS tiene los colores, eso se saca del archivo de configuración.

#### Tareas
* Andrea: Arreglar create_kar.py y escribir script para crear archivo CDSspositive.txt y CDSsnegative.txt
* Roberto
* Andrés

### Reunión 01/06
