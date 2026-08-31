# Instructions simples avant publication

Le dossier est déjà prêt pour GitHub avec le QMD, le script R, le rapport HTML, les trois CSV utiles, le petit RDS spatial et les huit cartes finales.

## Pour publier le dépôt

1. Décompresser l’archive `FEISTY_NWS_spatial_yield_analysis.zip`.
2. Sur GitHub, créer un dépôt vide, par exemple `FEISTY-NWS-spatial-yield-analysis`.
3. Ajouter tout le contenu du dossier décompressé au dépôt.
4. Vérifier que le fichier `README.md` apparaît bien sur la page d’accueil du dépôt.
5. Vérifier que les fichiers `nws_for_feisty_proj_year.RDS` et `all_FFMSY.RDS` ne figurent pas dans la liste publiée.

Les règles du fichier `.gitignore` empêchent normalement l’envoi accidentel de ces deux gros RDS.

## Pour rendre seulement le rapport

Tu n’as rien d’autre à copier. Ouvre `FEISTY_NWS.Rproj` dans RStudio, ouvre le QMD, puis clique sur **Render**. Les fichiers nécessaires au rapport sont déjà inclus.

## Pour relancer toute la simulation FEISTY

Copie localement les deux fichiers existants suivants dans le sous-dossier `data/` :

```text
nws_for_feisty_proj_year.RDS
all_FFMSY.RDS
```

Ils doivent donc se trouver exactement ici :

```text
FEISTY_NWS_spatial_yield_analysis/data/nws_for_feisty_proj_year.RDS
FEISTY_NWS_spatial_yield_analysis/data/all_FFMSY.RDS
```

Ne les publie pas sur GitHub. Le premier dépasse notamment la limite habituelle de taille d’un fichier GitHub.

Le petit fichier `NWS_full_grid_four_regions_cells.RDS` est déjà présent dans `outputs/rds/`, et `NWS_decadal_optimisation_best_approach.csv` est déjà présent dans `outputs/tables/`. Il ne faut ni les renommer ni les déplacer.

Enfin :

1. vérifie que le package R `FEISTY` est disponible dans ton environnement habituel ;
2. lance d’abord `NWS_FINAL_maps_FULL_22317_cells.R` avec `full_run <- FALSE` ;
3. si le test réussit, remplace cette valeur par `full_run <- TRUE` pour le calcul complet ;
4. rends de nouveau le QMD lorsque le calcul est terminé.

Le calcul complet est long, mais le script enregistre des points de reprise par blocs dans `outputs/rds/`.
