# Simulation de fluide SPH 2D

Simulation de fluide 2D par la méthode SPH (Smoothed Particle Hydrodynamics), réalisée en Processing (Java) pendant le master Informatique ISICG (Université de Limoges, 2024).

## Ce qui est implémenté

- SPH de base avec recherche de voisins par **grille**.
- **DFSPH** (Divergence-Free SPH) pour un fluide incompressible.
- **Turbulences** : confinement de vorticité (méthode classique), puis modèle micropolaire (rotation propre de chaque particule), d'après l'article [Turbulent Micropolar SPH Fluids with Foam](https://doi.org/10.1109/TVCG.2018.2832080) (Bender et al., 2019).
- Viscosité **XSPH**.
- Obstacles dans la scène, dont une hélice.

## Lancement

Processing exige que le dossier porte le nom du sketch : cloner le dépôt dans un dossier `fluid`, ouvrir `fluid.pde` dans [Processing](https://processing.org/) et lancer le sketch.

```bash
git clone https://github.com/JulesBaji/SPH-Fluid.git fluid
```
 Maintenir le clic de la souris ajoute des particules.

## Technologies

Processing (Java).
