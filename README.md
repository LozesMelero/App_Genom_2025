# App_Genom_2025
> *Répertoire qui contient les données sur les fichiers d'entrées pour quentiNemo, les sorties brutes ainsi que le code R pour les statistiues et graphiques associés*
> *Ce repertoire contient des données utiliser dans un cadre **pédagogique** de l'Université de Toulouse.*

---

## Description  des fichiers et dossier dans ce projet


### 1. Fichiers
- `fichiers.ini` : Fichiers utilisé pour paramétrer les simulations dans **quantiNemo** _(Neuenschwander, Michaud, & Goudet, 2019)_. 
- Les fichiers `_mean.txt`, `_var.txt`, `_legend.txt` et `_stats.txt` sont les fichiers de sorties, après avoir lancé les simulations.
- Les fichiers `.log` contiennent des informations sur les paramètres post simulations.

### 2. Dossiers
 - Le dossier `actions_plans` comporte les dossiers et fichiers de nos simulations en intégrant la passe à poisson sur le barrage A (`Remove_bar_A`) ou bien le barrage B (`Remove_bar_B`).
 - Le dossier `F_index` comporte les dossiers `Just_A`, `Just_B` et `Sans_barrage`. Ces dossiers ont été créés dans le cadre de notre scénario fictif où respectivement _seul le barrage A_, _seul le barrage B_ ont existé, et le dernier où _aucun des barrages_ n'a été construit. 

### 3. Script R 
- `App_Genom.R` script R utilisé pour les analyses statistiques et la création de graphiques se situe dans la dossier `Bin/Poisson`. Ce script à été consrtruit de sorte à pouvoir être utiliser en téléchargeant le présent répertoir GitHub.

 ![](./Logo_UT3.jpg =100x20)  ![](./truite_image.jpg) 
           
