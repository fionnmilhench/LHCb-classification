import sys
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

# Canvas settings
W, H = 9000, 3000
BG = 255

# Optionally use the seed
if len(sys.argv) > 1:
    seed = int(sys.argv[1])
    np.random.seed(seed)
else:
    seed = 2025

def raster_text(text, width_px=8):
    img = Image.new("L", (W, H), color=BG)
    draw = ImageDraw.Draw(img)
    size = int(H)
    try:
        font = ImageFont.truetype("/usr/share/X11/fonts/Type1/c0648bt_.pfb", size=size) # Bitstream Charter
    except (OSError, IOError):
        font = ImageFont.load_default()
    bbox = draw.textbbox((0, 0), text, font=font, anchor="la")
    tw = bbox[2] - bbox[0]
    th = bbox[3] - bbox[1]
    x = (W - tw) // 2
    y = -th // 4
    draw.text((x, y), text, fill=0, font=font)
    return np.array(img) < 128

def sample_mixture(mask, eps=0.05, n_total=2000000):
    rng = np.random.default_rng(seed)
    on_coords = np.argwhere(mask)
    k_text = int(round(eps * n_total))
    k_uni  = n_total - k_text
    perm = rng.permutation(on_coords.shape[0])
    chosen = on_coords[perm[:k_text]]
    jx = (chosen[:, 1] + rng.random(k_text)) / W
    jy = 1.0 - (chosen[:, 0] + rng.random(k_text)) / H
    pts_text = np.column_stack([jx, jy])
    pts_uni  = rng.random((k_uni, 2))
    pts = np.vstack([pts_text, pts_uni])
    rng.shuffle(pts)
    return pts

def show_density(pts, title, bins=400):
    H2, xe, ye = np.histogram2d(pts[:,0], pts[:,1], bins=bins, range=[[0,1],[0,1]])
    H2 = H2.T
    plt.figure()
    plt.xlabel(r'$\zeta_1$')
    plt.ylabel(r'$\zeta_2$')
    plt.imshow(H2, origin='lower', extent=[0,1,0,1], aspect='equal', cmap='BuGn')
    plt.colorbar(label='counts per bin')
    plt.title(title)

# Generating the text-encoded distribution
mask_lhcb = raster_text("LHCb")
pts_lhcb = sample_mixture(mask_lhcb, eps=0.07, n_total=2000000)
show_density(pts_lhcb, "Binned text ensemble")
np.savetxt("datasets/encoded_text_distribution_lhcb_ensemble.txt", pts_lhcb)
plt.savefig("plots/encoded_text_distribution_lhcb_ensemble.png", dpi=300)

# One sample dataset of text-encoded points
pts_lhcb_sample = pts_lhcb[:20000]
show_density(pts_lhcb_sample, "Binned sample text distribution", bins=200)
np.savetxt("datasets/encoded_text_distribution_lhcb_single.txt", pts_lhcb_sample)
plt.savefig("plots/encoded_text_distribution_lhcb_single.png", dpi=300)

# Generating the uniform random distribution
pts_random = np.random.default_rng(seed).random((2000000, 2))
show_density(pts_random, "Binned random ensemble")
np.savetxt("datasets/distribution_random.txt", pts_random)
plt.savefig("plots/distribution_random.png", dpi=300)

# One sample dataset of random points
pts_random_sample = pts_random[:20000]
show_density(pts_random_sample, "Binned sample random distribution", bins=200)
np.savetxt("datasets/distribution_random_single.txt", pts_random_sample)
plt.savefig("plots/distribution_random_single.png", dpi=300)