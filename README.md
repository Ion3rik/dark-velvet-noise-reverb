# DARK-VELVET-NOISE REVERB
Code for the dark-velvet-noise artificial reverberation algorithm [1].

# Current Features
- Offline Matlab implementation, including the full pipeline from fitting the EDVN model parameters to a target late-reverberation impulse response to generating the model IR. See ```example_IrModel.m```.
- Binaural extension that allows fitting the model to binaural late-reverberation with matched interaural coherence [2]. See ```example_BrirModel.m```. 
- Matlab Demo app, providing example parametrization for modifying the revebreration. See ```dvn_reverb.mlapp```.

# Upcoming Features
- Python version of the core Matlab code.
- Real-time plugin implemented with Juce and C++

# References
1. Jon Fagerström, Sebastian J. Schlecht, Vesa Välimäki. Non-Exponential Reverberation Modeling Using Dark Velvet Noise. Journal of the Audio Engineering Society, Vol. 72, 6, pp. 370–382, June 2024. https://aes2.org/publications/elibrary-page/?id=22636
2. Jon Fagerström, Nils Meyer-Kahlen, Sebastian J. Schlecht, Vesa Välimäki. Binaural Dark-Velvet-Noise Reverberator. In Proceedings of the International Conference on Digital Audio Effects (DAFx), pp. 246–253, September 2024. https://www.dafx.de/paper-archive/2024/papers/DAFx24_paper_63.pdf
