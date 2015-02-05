#!/bin/bash
if (($#==2))			#verification du nb d'argument
then
	if [ -d "./archives_figures/$1" ] 		#verification que le fichier n'existe pas
	then
		echo "ERREUR: le fichier existe deja!"
		echo "voulez-vous le remplacer?(O/N)"
		read choice
		if [[ $choice == "O" ]]
		then   
			echo "suppression du dossier de figure nommé $1 ."
			rm -R ./archives_figures/$1
			echo "Création du nouveau repertoire de sauvegarde...."
		else
			echo "Arret du programme"
			exit 1
		fi
	else
	echo "Création du repertoire de sauvegarde...."
	fi
		touch log.txt   													   			
		touch $PPID.m                                                                             #creation d'un script .m temporaire
		cat main.m > $PPID.m                                                                      #reprise du code complet
		echo "for (i=1:$2) saveas(i,['figure'  num2str(i) '.png'],'png'); end" >> $PPID.m         #ajout des lignes de print en png

#=======================================================================================================================
#Sauvegarde des coefficients de reflexions et transmission en txt
#        echo "coll1 = ones(N,1); coll2 = ones(N,1); coll3 = ones(N,1); coll4 = ones(N,1);">>$PPID.m
#		echo "for (i =1:N) coll1(i) = real(R)(1,1,i); coll2(i) = imag(R)(1,1,i); coll3(i) = real(T)(1,1,i); coll4(i) = imag(T)(1,1,i); end" >> $PPID.m
#		echo "T_sauvegarde = [f' coll3 coll4]; R_sauvegarde = [f' coll1 coll2]; " >> $PPID.m
#		echo "save('R.txt', 'R_sauvegarde'); save('T.txt', 'T_sauvegarde');" >>$PPID.m 
#======================================================================================================================
	
		octave --silent --no-window-system $PPID.m { >> log.txt }   	#éxecution du script modifié et récupération du résultat

		mkdir $1	#creation du dossier de reception

		mv log.txt $1			#remplissage du dossier de reception
		for i in $(seq $2); do
			mv figure$i.png $1
		done

		mv $1 archives_figures			#rangement du dossier de reception
		rm $PPID.m 		

#======================================================================================================
#		mv R.txt $1
#		mv T.txt $1															
#==============================================================================================                              						                 
else
  echo "ERREUR: Il y a un probleme d'argument => mettre un mot et le nb de figure"	#erreur pour argument différent de 2
fi

echo "Figures archivées"
exit 0
