# GenomeVisualizer

Genome visualizer software project

### Reunión 06/01

* GBK a fna, GBK a faa, GBK a gff.
* Mirar proyecto [SnakemakeCircos](https://bitbucket.org/mroachawri/snakemakecircos/src/master/). Ver cómo se podría adaptar.
* Ver cómo mejorar la sintaxis.
* Ver qué cosas se puede agregar (Cog).
* Hacer grupo de Telegram

### Reunión 04/05

* Generalmente los archivos están en GBK y usar este archivo es lo mejor ya que te da toda la información que otros archivos te dan por partes (fna, faa y gff).
* Nos importa graficar contenido GC, simetría GC, secuencias codificantes que encontró (verde una hebra y azul la otra) y tamaño del genoma (cada fragmento representa un contig).
* Opcionalmente si el usuario quiere hacer la clasificación COG (que aveces viene en GBK), se necesitan todas las secuencias de nucleótidos, por lo que ahí sale más fácil que el usuario ingrese el archivo fna. Para obtener la anotación se puede correr un [DFAST](https://dfast.ddbj.nig.ac.jp/) con el archivo fna. Esto produce un nuevo GBK que sí tiene información sobre COG categories.
* chr.kar es un archivo input para CIRCOS. Este archivo se crea tomando el valor indicado en el campo SOURCE del archivo gbk. La primera línea contiene el tamaño total del genoma. "band" simboliza cada uno de los contigs. Un ejemplo de SOURCE para dos contigs consecutivos: SOURCE_1 "1...N" y SOURCE_2 "1...M". Para ponerlo en chr.kar se debe tomar "0...N" y "N+1...N+M".
* Andrea: Tener esto listo para la próxima semana.
* Roberto: Revisar lo que hace DFAST.
* Se podría también graficar los repeat_regions, tRNA, rRNA. Además de destacar ciertas regiones que el usuario puede elegir.
* Vicente se va a unir al proyecto.

### Reunión 18/05


