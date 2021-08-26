# diretorio com arquivo de entrada
dirExp="$(pwd)"/${1:-${dirPadrao}} 

#### DEFINICAO DO EXECUTAVEL ####
LOCAL=$(pwd)
NOMEEXECUTAVEL=run
DIRBIN="${LOCAL}"
#rm -f ${LOCAL}/output.out

EXECUTAVEL=${DIRBIN}/${NOMEEXECUTAVEL}

#### definicao do comando a ser executado
maquina=$(hostname)
comando="( cd ${dirExp}; time  ${EXECUTAVEL})"

if [ -e ${EXECUTAVEL} ] 
then
  printf "\n diretorio do experimento.: %s\n" ${dirExp}  
  printf "\n nome do executavel.......: %s\n" ${EXECUTAVEL}
  printf "\n comando .................: %s\n" "${comando}"
  eval ${comando} 
#  eval ${comando} > output.out
#  eval ${comando}  |tee  ${arqTela}
else
  printf "\n EXECUTAVEL NAO ENCONTRADO \n"
  printf "\n comando .................: %s\n" "${comando}"
fi
