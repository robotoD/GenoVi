import requests
import xmltodict

# Parse assembly xml to get all human chromosomes
with open("./human.xml", "rb") as xml:
    human = xmltodict.parse(xml)
accessions = []
human = human['ROOT']['ASSEMBLY']['CHROMOSOMES']['CHROMOSOME']
for i in range(len(human)):
    accessions.append([human[i]['NAME'], human[i]['@accession']])


# Use requests to download all the chromosomes


def ena_download(list):
    url = "https://www.ebi.ac.uk/ena/data/view/{}&display=fasta".format(list[1])
    r = requests.get(url, stream=True)
    with open("{}.fasta".format(list[1]), "wb") as file:
        for chunk in r.iter_content(chunk_size=1024 * 512):
            if chunk:
                file.write(chunk)


if __name__ == "__main__":
    for i in accessions:
        ena_download(i)