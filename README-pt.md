parsnv é uma aplicação para executar em paralelo um pipeline de chamada de
variantes usando samtools, bcftools, e vcftools, visando baixo uso de memória.

# Instalação

A forma mais fácil de instalar e usar o parsnv é utilizando [esta imagem
docker](https://hub.docker.com/repository/docker/afarah1/ubuntu-samtools), que
fornece um ambiente pronto para execução do programa. Nela, basta baixar o
arquivo parsnv desse repositório e torná-lo executável:

```
wget https://github.com/afarah1/parsnv/blob/main/parsnv && chmod +x parsnv
```

Caso você não queira ou não possa usar docker, certifique-se de que as
dependências abaixo tenham sido cumpridas. Em caso de problemas consulte a
seção "Erros comuns".

## Dependências

Você deve possuir GNU Parallel, samtools, bcftools, e vcfutils no seu PATH.
Consulte o manual dessas ferramentas para realizar sua instalação.

## Docker

Para obter a imagem basta usar [pull](https://docs.docker.com/engine/reference/commandline/pull/). Para não precisar copiar os dados para dentro da imagem docker, utilize [volumes](https://docs.docker.com/storage/volumes/) ou [bind mounts](https://docs.docker.com/storage/bind-mounts/). Por exemplo, se seus dados estão em `/DATA`, você pode obter a imagem e usar bind mounts da seguinte forma:

```
docker pull afarah1/ubuntu-samtools
docker run -it --mount type=bind,source=/DATA,target=/DATA afarah1/ubuntu-samtools
```

# Uso

`parsnv [opções] lista`

`lista` é um arquivo de texto contendo o caminho completo para os arquivos de
sequências (CRAMs) sobre os quais deseja-se chamar as variantes, bem como
arquivos do genoma de referência. Um exemplo é fornecido abaixo, bem como uma
descrição completa do conteúdo esperado do arquivo. A saída da aplicação é
escrita para um arquivo chamado `final.vcf` contendo as variantes anotadas.

## Conteúdos da lista

### Arquivos CRAM 

São esperados arquivos com extensão `.cram` - as amostras de sequências genéticas
alinhadas. Somente CRAM é suportado, mas o código é fácil de adaptar para
suportar BAM ou SAM.

Se você tiver arquivos de índice no mesmo diretório e com o mesmo nome de
arquivo `.cram`, mas com `.crai` ao final, eles serão reconhecidos e as
sequências não serão re-indexadas, poupando tempo (não muito). Não é necessário
listar os `.crai`, basta estarem no mesmo diretório dos `.cram`.

### Arquivo FASTA 

É esperado um único arquivo FASTA - o genoma de referência ao qual as
sequências estão alinhadas. A extensão deve ser `.fa` ou `.fasta`, podendo estar
comprimido como `.gz`.

### Arquivo VCF 

É esperado um arquivo contendo as variantes do genoma de referência. A extensão
deve ser `.vcf.gz` (deve estar comprimido).

### Índice VCF

É esperado o arquivo de índice do VCF supracitado. A extensão deve ser
`.vcf.gz.tbi` ou `.vcf.gz.csi`.

Se você não fornecer esse arquivo, o programa irá tentar baixá-lo do NCBI (o
arquivo possui poucos MBs).

## Exemplo

Abaixo está um exemplo de lista para chamada de variantes nas amostras
HGDP00525, HGDP00519, HGDP00711, HGDP00712, HGDP00714, HGDP00715, HGDP00719,
e HGDP00777, alinhadas ao GRCh38.

```
/data/HGDP00525.alt_bwamem_GRCh38DH.20181023.French.cram
/data/HGDP00519.alt_bwamem_GRCh38DH.20181023.French.cram
/data/HGDP00711.alt_bwamem_GRCh38DH.20181023.Cambodian.cram
/data/HGDP00712.alt_bwamem_GRCh38DH.20181023.Cambodian.cram
/data/HGDP00714.alt_bwamem_GRCh38DH.20181023.Cambodian.cram
/data/HGDP00715.alt_bwamem_GRCh38DH.20181023.Cambodian.cram
/data/HGDP00719.alt_bwamem_GRCh38DH.20181023.Cambodian.cram
/data/HGDP00777.alt_bwamem_GRCh38DH.20181023.Han.cram
/data/GRCh38_full_analysis_set_plus_decoy_hla.fa
/data/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
/data/00-All.vcf.gz
/data/00-All.vcf.gz.tbi
```

# Opções

O programa permite as seguintes opções:

```
  -i, --indels
    Não remove indels do VCF

  -f, --filter arquivo
    Lê de um arquivo de texto (passado como parâmetro) um filtro customizado
    para o vcftools. O filtro padrão é --max-missing 0.5 --minDP 5
    --min-alleles 2 --max-alleles 2 --minQ 20.

  -h, --help
    Exibe um texto de ajuda (em inglês).
```

# Erros comuns

Esta seção descreve erros comuns e como resolvê-los.

## Erros do SAMtools 

### Índice corrompido

```
[main_samview] retrieval of region \"%s\" failed due to truncated file or corrupt BAM index file
```

Isso significa que uma leitura do samtools falhou, provavelmente devido a um
arquivo de índice corrompido. Isso é um erro do samtools, não do parsnv. Tente
remover todos arquivos `.crai` e rodar o programa de novo.

### CRAM corrompido

```
samtools index: failed to create index for "sample.cram": No such file or directory
```

Isso significa que o samtools falhou ao tentar indexar um CRAM. Isso é um erro
do samtools, não do parsnv. A mensagem de erro emitida pelo samtools sugere que
ele não encontrou o arquivo, mas na verdade o samtools emite ela quando a
indexação falhou por um motivo desconhecido, não necessariamente porque ele não
encontrou o arquivo. O parsnv verifica que todos arquivos da lista estão
acessíveis, então a não ser que você tenha movido os arquivos ou mudado suas
permissões após já ter iniciado o parsnv, provavelmente o problema é que o CRAM
em questão está corrompido. Isso geralmente é causado por um download
incompleto. Verifique o tamanho do seu CRAM contra o tamanho do CRAM no FTP de
onde você obteu ele, se os tamanhos não são iguais seu download ficou
incompleto e você precisará baixar o CRAM de novo. Se você está obtendo os
CRAMs de uma fonte que fornece um checksum compare o checksum e não o tamanho
dos arquivos.

## Erros do parsnv 

### Comando não encontrado

```
parallel/samtools/bcftools/vcftools: command not found
```

Isso significa que você não possui a ferramenta em questão ou ela não está no
seu PATH. Vide a seção Instalação.

### Arquivo não encontrado

```
Could not find FILE on cram list
```

A lista está incompleta. Vide a seção "Uso".

### Arquivo inacessível

```
Could not access file FILE
```

O arquivo fornecido não está acessível. Verifique o seguinte:

1. O caminho está correto (teste com `ls`)
2. O caminho está completo (não é um caminho relativo ao diretório atual)
3. O arquivo está acessível (verifique as permissões dele)

# Pipeline

O pipeline será descrito abaixo, por etapas. As saídas de cada comando não são
transmitidas ao próximo via pipe, são escritas no disco e os arquivos
intermediários são deletados quando não mais necessários. Isso é proposital,
mas pode deixar a execução mais lenta do que se pipes fossem usados. O
benefício é um uso baixo de memória.

As opções `filter` e `indels` permitem uma parametrização rudimentar. Além
disso, o pipeline está escrito de forma a atender as necessidades de um usuário
em particular. Você pode precisar ajustar o código para suas necessidades
particulares, caso elas não sejam atendidas através das opções existentes.

## Etapa 1

Indexa os CRAMs e separa-os em 22 regiões:

```
FOR CRAM_J IN LIST DO PARALLEL
  INDEX CRAM_J
  FOR REGION_I IN 1,22 DO
    VIEW CRAM_J REGION_I >CJ_RI
```

## Etapa 2

Aglutina os CRAMs por região, executa a chamada de variantes para cada região:

```
FOR REGION_I IN 1,22 DO PARALLEL
  MERGE (CJ_RI FOR J IN LIST)>MERGED_I
  INDEX MERGED_I
  MPILEUP
  CALL
  BCF VIEW
  VCF REMOVE INDELS
  VCF FILTER
```

## Etapa 3

Aglutina sobre todas regiões, anota, e filtra as variantes:

```
VCF CONCAT
BCF VIEW
BCF INDEX
BCF ANNOTATE
BCF VIEW
VCF FILTER
```
