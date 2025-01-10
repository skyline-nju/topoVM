import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image as mpimg


if __name__ == "__main__":
    fig, axes = plt.subplots(1, 3, figsize=(8, 3), sharex=True, sharey=True)

    snaps = ["fig/snap/alpha_L400_D0.03.png",
        "fig/snap/alpha_L400_D0.20.png",
        "fig/snap/alpha_L400_D0.31.png"
        ]
    
    titles = [
          r"(a) $D=0.03$",
          r"(b) $D=0.20$",
          r"(c) $D=0.31$",
    ]
    
    for i, ax in enumerate(axes):
            ax.set_title(titles[i], fontsize="xx-large")
            im = mpimg.imread(snaps[i])
            ax.imshow(im)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal")

    dx = 0.25
    ax_cb = axes[1].inset_axes([0, 0, dx, dx])
    ax_cb.set_xticklabels([])
    ax_cb.set_yticklabels([])
    ax_cb.set_xticks([])
    ax_cb.set_yticks([])
    im = mpimg.imread("fig/circle2.png")
    ax_cb.imshow(im)
    # ax_cb.set_title(r"$\theta_i$", fontsize=fs)
    ax_cb.set_title(r"$\theta_i$", fontsize="xx-large")
    fig.tight_layout(rect=[-0.015, -0.03, 1.015, 0.97], w_pad=0.5)
    plt.show()
    # plt.savefig("fig/FIG5.pdf", dpi=150)
    plt.close()