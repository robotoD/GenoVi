# GenomeVisualizer

Genome visualizer software project


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
* Roberto: Revisar lo que hace DFAST. Ver cómo se anotan los COGs. Preguntar a Vicente si se quiere unir al proyecto.
* Andrés: Averiguar qué información se debe extraer del GBK para hacer los archivos input para colorear secciones custom.


### Reunión 18/05

* Prueba de primera versión de create_kar.py
* De create_kar modificar: El archivo .kar final debe quedar en el mismo orden de como aparecen los contigs en el .gbk. El orden de las bandas se cuenta desde 0 de nuevo. Agregar prefijos de argumentos, verificadores y comentarios a create_kar.
* Siguiente archivo a crear es el que define las flechas de las secuencias codificantes (CDSs_positive.txt y CDSs_negative.txt). En el archivo gbk, en cada contig aparece el feature CDS, el location indica donde inicia y termina la secuencia. Si aparece complement(2..444) significa que la secuencia es negativa. Si aparece 2..444 la secuencia es positiva.
* A la locación de todos los CDS de un contig se les debe sumar el tamaño total del contig anterior --El orden de los contigs da igual--. Por ejemplo, si dice el contig anterior dice "source 1..10000" y el primer CDS del contig siguiente tiene location "4..10", debe quedar como "10004..10010".
* La última columna del CDS tiene los colores, pero esta informacion no es necesaria. Esto se grafica desde el archivo de configuración.
* Pipeline para la anotacion de COGs y obtencion de un una tsv con la informacion a graficar: (https://github.com/transcript/COG). Requiere descargar la base de datos COGs del NCBI, Instalar DIAMOND, realizar la base de datos de DIAMOND, hacer la búsqueda con DIAMOND en el .faa. Entrega una tabla con 3 columnas, la tercera es el código COG). Cada uno de los 25 grupos COGs (e.g. A,B,C,D,E...) debe ser ingresado como un archivo aparte (COG_K.txt) indicando el inicio y final de cada CDS clasificado en cada COG (indicado en el archivo anterior). 
* Vicente esta abierto a participar.

#### Tareas
* Andrea: Arreglar create_kar.py y escribir script para crear archivo CDSspositive.txt y CDSsnegative.txt
* Roberto: Corroborar pipeline de anotacion COGs para integrarlo al programa. Buscar otros analisis que agregar al visualizador.
* Andrés: Verificar el pipeline para obtencion de datos de GC-skew y GC-content desde un archivo ###.fna para crear el script de python y automatizar la creacion de archivos.
* Vicente: Script de python obtener desde un ###.gbk, un archivo ###.faa (para análisis de COGs) y un archivo ###.fna (para cálculo de GC-Skew y GC-content).  

### Reunión 01/06
* Vicente Saona se une al proyecto (Wohooo!).
* Archivos CDSpositive y CDSnegative tienen que estar separados por tabs. La extensión da lo mismo. Además, la location de los contigs en los CDS debe llevar la suma total de todos los contigs anteriores (como en el kar). La primera columna de los ODS es la misma primera columna del kar.
* Hay que editar el gbk que se va a transformar a fna

#### Tareas
* Andrea: 
    1. Arreglar creación de CDS. 
    2. Apoyar a Roberto con lo de los COGs.
* Roberto: Continuar con el pipeline de anotación COGs.
* Andrés: Continuar viendo que se hace con los archivos de salida, juntarlos, etc, para obtener los archivos de configuración de Circos. 
* Vicente: Modificar los scripts que encontró Andrés para incorporar los números que le faltan para el inicio y el fin, con tal de que el output quede listo para Circos.

### Reunión 22/06
* Andrés hizo nuevas figuras usando los scripts. Le facilitaron la pega pero aún falta agregar cosas.
* Idea: crear paletas predefinidas que el usuario pueda elegir para usar en sus visualizaciones.

#### Tareas
* Andrea: 
    1. Cambiar "chr1" en vez de "chr01". 
    2. Crear archivos trnas, sería lo mismo que los CDSs, pero buscando "trna"s. 
    3. Apoyar a Roberto con lo de los COGs.
* Roberto: Continuar con el pipeline de anotación COGs.
* Andrés: 
    1. Pensar qué otro análisis podemos hacer, tal vez agregar una funcionalidad para graficar genomas cerrados. 
    2. ¿Cómo agregarle texto a las imágenes generadas por circos?. 
    3. Si agregas un chr2, ¿lo graficaría igual?. 
    4. (A futuro) Tomar notas de las cosas que se tienen que hacer manuales cuando se usan los scripts para crear visualizaciones.
* Vicente: 
    1. En GC content la primera columna debe decir "chr1" en todas las filas. 
    2. Eliminar la separación entre contigs en el archivo fna (borrando los ">n"). 
    3. Buscar mínimo y máximo de los GC content para definir el mínimo y máximo en la visualización.

### Reunión 20/07
* Andrés trabajó en poner texto en la imagen, ver qué pasaba si agrega un segundo cromosoma al visualizador y jugar con los valores de GC content. Roberto sugiere cambiarlo para que tome varios cromosomas y haga la separación. Habría que modificar la creación de kars y CDSs.
* Vicente hizo sus cambios, faltaría solamente agregar más de un cromosomas.
* Pensar en cómo agregar paletas de colores para que el usuario pueda elegir.

#### Tareas
* Andrea: 
    1. Crear archivos trnas, sería lo mismo que los CDSs, pero buscando "trna"s. 
    2. Creación de KAR y CDSs debe aceptar más de un cromosoma.
    3. Apoyar a Roberto con lo de los COGs.
* Roberto: Continuar con el pipeline de anotación COGs.
* Andrés:
    1. Identificar cuáles son los parámetros variables del archivo config de Circos.
    2. Poner positivo y negativo el GC-skew.
* Vicente:
    1. Considerar más de un cromosoma en el GC analysis.
    2. Ver si se le pueden poner títulos al centro o arriba a los SVG.
    3. Ver cómo generar múltiples círculos individualmente y ver cómo se pueden juntar pero con tamaños a escala de la cantidad de kb.

### Reunión 03/08
* Roberto encontró un programa basado en machine learning llamado DeepNOG, que predice las categorías de los COGs.
* Andrés logró hacer que el GC-skew sea de distinto color para positivo y negativo. Debe descubrir cómo cambiar los colores positivo y negativo porque ahora están fijos en celeste y azul. Tal vez se podría hacer que las ventanas sean más anchas. Identificó los parámetros del archivo conf que puede variar de genoma en genoma. 
* Vicente logró hacer la generación del archivo GC content soporte más de un cromosoma. Además hizo que a las imágenes se les pueda agregar un título. También logró hacer que las imágenes puedan juntarse con diferente tamaño.
* Una mejora a futuro para la creación de los GC contents es que el script pueda recibir archivos de contigs independientes.

#### Tareas
* Andrea: 
    1. Creación de archivos debe soportar más de un cromosoma.
    2. Dividir archivo final de COGs en positivo-negativo y por categoría (p. ejemplo A_positivo y A_negativo). 
* Vicente:
    1. Hacer script en Python que una todo. 
    2. Hacer que las imágenes juntadas queden proporcional al tamaño de los genomas.
* Andrés: 
    1. Averiguar cómo cambiar los colores de los GC-skew positivo y negativo.
    2. Poner en un word y marcar los parámetros del conf que varían.

### Reunión 31/08
* Andrés probó el script completo y se percató de cosas que modificar el los archivos de salida. Entre ellos el kar.file (corregir el inicio y término de cada "chr" y "band", dejar los "chr" siempre en "black" y las bandas intercalarlas en "white" y "black") y los cds_pos cds_neg (rectificar los locations).  
* Además 
* Vicente editó el archivo genbank2faa.py para transformar archivos gbk/gbff que tenian problemas. Uno de los problemas detectados era en archivo gbk/gbff de NCBI que incluye anotacion de pseudogenes, estos eran identificados como CDS pero sin feature de translation (al ser pseudogen no se anota). La edicin hace que en estos casos el script se salta el CDS. 
* Andrea integro DeepNOG al analisis para realizar la clasificacion de COGs a las CDS extraidas del gbk/gbff file. Esto genera 46 archivos de categories 23 de la hebra positiva y 23 de la hebra negativa. Estos deben ser incorporados al highlights.conf para ser graficados por Circos.

#### Tareas
* Andrea: 
    * Editar el create_raw.py para incluir las modificaciones: 
    1. Crear siempre el archivo kar file donde cada contig es un "chr" y "band" independiente 
    2. Editar el kar file: corregir el inicio y término (locations) de cada "chr" y "band"; dejar los "chr" siempre en "black" y las bandas intercalarlas en "white" y "black"). 
    3. Editar la creacion de los archivos del cds_pos y cds_neg (rectificar las locations de cada CDS)
    4. --complete_genome debera separar el archivo original gbk/gbff (1 file) en archivos "n" dependiendo de la cantidad de contigs que posea y realizar el pipeline completo a cada gbk/gbff creado por separado. Luego realizar el ajuste de tamaño realizado por Vicente (mergeImages.py)
* Vicente:
    1. Editar genbank2faa.py para que al extraer cada secuencia guarde el locus_tag y el filename del gbk/gbff como identificador de cada secuencia
    2. Editar el archivo svg para incluir cuadrados con los colores correspondientes a cada clasificacion (CDS, GC-content, GC-skew, COGs) con su leyenda.
* Andrés: 
    1. Agregar en el archivo highlight.conf la posicion de los COGs categories
* Roberto:
    1. Identificar los RGB para cada uno de los COGs anotados por DeepNOG (colores por macrocategorias)

### Reunión 21/09
* Andrés tuvo problemas con tu PC, sufrio mucho. No pudo utilizar los scripts, en general no identificaba el 1er contigo. Era un problema con algunos gbk que vienen de la anotacion de prokka, donde el ID de cada contigo puede variar y genera conflicto.
* Vicente editó el genomeVisualizer.py para integrar todas las partes a excepcion del análisis de deepnog. Funciona y genera la visualizacion automatica utilizando circos. Esto sin incorporar el comando (-gc)
* Andrea modifico el create_kar.py, ahora realiza de forma correcta el archivo ####bands.kar y los cds_pos.txt y cds_neg.txt (esto al usarlo sin el comando -gc).
* Roberto realizo una escala de color en RGB para tener colores unicos para los COGs a graficar, fueron separados en clases mayores las cuales comparten colores similares.


#### Tareas
* Tod@s: Revisar los colores para ver modificaciones en la escala de colores de COGs.
* Andrea: 
    * Editar el create_raw.py para incluir las modificaciones: 
    1. Modificar la identificacion del nombre de cada contig (ID_contig) para no tener problemas con algunos archivos (especialmente PROKKA) 
    2. Revisar que el create_kar.py para que no elimine la primera linea de cada gbk luego de correr.
    3. Crear los archivos cds_pos.txt y cds_neg.txt al usarlo utilizando el comando -gc --get_categories
    4. Integrar la creacion de los archivos rrna (features).
* Vicente:
    1. Incorporar el pipeline de deepnog en el genomeVisualizer.py para realizar todo inmediatamente
* Andrés: 
    1. Editar el archivo highlight.conf para incorporar: radius length de cada feature (tRNA, rRNA, CDS_pos, CDS_neg, COGs_pos, COGs_neg), pathfile de cada archivo y colores a utilizar.
* Roberto:
    1. Identificar un valor de corte para la integracion de DeepNOG (hay identificacion de proteinas con un valor de confidence muy bajo)
