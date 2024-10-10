# Fichier README

## Architecture du Code

Le code est organisé en trois dossiers :

- **Src**: Contient le code source.
- **Results**: Contient les résultats après exécution.
- **Data**: Contient les fichiers de données (.toml).

## Notice du Code

Pour lancer une simulation, suivez les étapes suivantes :

1. **Préparation du fichier de données**:
   - Ouvrir un fichier 'data_cas_test_x.toml' dans le dossier 'Data'.
   - Choisissez le nombre de points 'Nx' et 'Ny', la 'CFL', le schéma en espace 'space_scheme' (1 ou 2), le schéma en temps 'time_scheme' (1 ou 2) et le pas de temps 'dt'. Notez que seules les combinaisons de schémas suivantes sont possibles : (1, 1) et (2, 2).

2. **Nettoyage du code source** :
   - Ouvrez un terminal dans le dossier 'Src'.
   - Tapez la commande 'make clean'.

3. **Compilation du code** :
   - Toujours dans le même terminal, tapez la commande 'make'.

4. **Exécution du code** :
   - Dans le même terminal, exécutez le code en passant en argument un fichier '.toml'.
   - Exemple: 'mpirun -n 4 ./run ../Data/data_cas_test_1.toml'.

5. **Traitement des résultats** :
   - Ouvrez un terminal dans le dossier 'Results'.
   - Fusionnez les fichiers de sortie.

6. **Visualisation** :
   - Toujours dans le même terminal, tapez la commande 'gnuplot visu.gnu'.

