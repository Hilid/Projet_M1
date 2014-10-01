#!/bin/bash
if (($#==2))																		#verification du nb d'argument
then
	if [ -d "./archives_figures/$1" ] 												#verification que le fichier n'existe pas
	then
		echo "ERREUR: le fichier existe deja!"
	else
		cp config.txt log.txt   													   #sauvegarde de la configuration du réseau
		touch $PPID.m                                             					     #creation d'un script .m temporaire
		cat main.m > $PPID.m                                   								        #reprise du code complet
		echo "for (i=1:$2) saveas(i,['figure'  num2str(i) '.png'],'png'); end" >> $PPID.m            #ajout des lignes de print
		octave --silent --no-window-system $PPID.m { >> log.txt }   					#éxecution du script modifié et récupération du résultat

		mkdir $1																		#creation du dossier de reception

		mv log.txt $1																	#remplissage du dossier de reception
		for i in $(seq $2); do
			mv figure$i.png $1
		done

		mv $1 archives_figures															#rangement du dossier de reception
		rm $PPID.m 																		#suppression du fichier temporaire
	fi                              						                 
else
  echo "ERREUR: Il y a un probleme d'argument => mettre un mot et le nb de figure"			#erreur pour argument différent de 2
fi

exit 0
