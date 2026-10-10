"""Noise reduction in astro-images: mean vs median vs Gaussian filters, with SNR.

SNR = mean star peak / standard deviation of the star-free background.
By default a synthetic star field with known stars is used, so peak preservation can be
measured exactly. Use --image to run the same filters on your own photograph.

    python spatial_filtering.py
    python spatial_filtering.py --image "Omega Centauri.JPG"
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage as ndi


def synthetic_field(size=256, n_stars=40, seed=0):
    """Gaussian stars + Gaussian read noise + hot pixels (impulse noise)."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[:size, :size]
    clean = np.zeros((size, size))
    centres = rng.integers(8, size - 8, size=(n_stars, 2))
    for (cy, cx), amp in zip(centres, rng.uniform(0.4, 1.0, n_stars)):
        clean += amp * np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * 1.5**2))
    noisy = clean + rng.normal(0, 0.08, clean.shape)
    hot = rng.random(clean.shape) < 0.01
    noisy[hot] = rng.uniform(0.5, 1.0, hot.sum())
    return clean, np.clip(noisy, 0, 1)


def masks_from(reference):
    """Signal = brightest pixels; background = pixels far from any bright region."""
    smooth = ndi.median_filter(reference, size=5)
    signal = smooth > np.percentile(smooth, 99)
    background = ~ndi.binary_dilation(smooth > np.percentile(smooth, 90), iterations=6)
    return signal, background


def snr(img, signal_mask, background_mask):
    return img[signal_mask].mean() / img[background_mask].std()


def apply_filters(img):
    return {
        "original": img,
        "mean 3x3": ndi.uniform_filter(img, size=3),
        "median 3x3": ndi.median_filter(img, size=3),
        "Gaussian sigma=1": ndi.gaussian_filter(img, sigma=1.0),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--image")
    ap.add_argument("--out", default="spatial_filtering.png")
    a = ap.parse_args()

    if a.image:
        from PIL import Image
        noisy = np.asarray(Image.open(a.image).convert("L"), dtype=float) / 255.0
        clean = None
    else:
        clean, noisy = synthetic_field()

    signal_mask, background_mask = masks_from(clean if clean is not None else noisy)
    results = apply_filters(noisy)

    print(f"{'filter':<18} {'SNR':>7} {'bg std':>8} {'star peak':>10}")
    for name, img in results.items():
        print(f"{name:<18} {snr(img, signal_mask, background_mask):7.2f} "
              f"{img[background_mask].std():8.4f} {img[signal_mask].mean():10.3f}")

    fig, axes = plt.subplots(1, 4, figsize=(16, 4.2))
    for ax, (name, img) in zip(axes, results.items()):
        ax.imshow(img, cmap="gray", vmin=0, vmax=1)
        ax.set_title(f"{name}\nSNR = {snr(img, signal_mask, background_mask):.1f}")
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)


if __name__ == "__main__":
    main()
