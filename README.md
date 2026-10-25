# Simulation de fluide SPH 2D

Simulation de fluide 2D par la méthode SPH (Smoothed Particle Hydrodynamics), réalisée en Processing (Java) pendant le master Informatique ISICG (Université de Limoges, 2024).

## Ce qui est implémenté

- SPH de base avec recherche de voisins par **grille**.
- **DFSPH** (Divergence-Free SPH) pour un fluide incompressible, d'après l'article de recherche.
- **Confinement de vorticité** (modèle classique) et **modèle micropolaire**.
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
