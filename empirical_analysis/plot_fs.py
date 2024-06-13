# micromamba activate low-coverage

import os
import dadi
from matplotlib import pyplot as plt

# import data
angsd_fs_dir='YRI_ANGSD/archive_GL_GATK/angsd_results'
all_data_fs, all_angsd_fs = [], []

for i, depth in enumerate(('3x', '5x','10x', '30x')):
    data_fs = dadi.Spectrum.from_file(f'YRI/inference_two_epoch/{depth}_fs_subsampled_32')
    all_data_fs.append(data_fs)

    fn=f"{depth}.sfs.mle.txt"
    # load ANGSD fs
    with open(os.path.join(angsd_fs_dir, fn)) as fh:
        fs_txt = fh.readline().strip()
    fs = [float(entry) for entry in fs_txt.split()]
    angsd_fs = dadi.Spectrum(fs, mask_corners=True)
    
    all_angsd_fs.append(angsd_fs)

# plotting
fig = plt.figure(figsize=(14, 5))
plt.rcParams.update({'font.size': 12})

ax1 = fig.add_subplot(121)
ax1.set_yscale("log")
plt.tick_params('both', length=7, which='major')
ax2 = fig.add_subplot(122, sharey=ax1)
plt.tick_params('both', length=7, which='major')

ax1.set_title("GATK Allele Frequency Spectra", fontsize=12)
ax1.plot(all_data_fs[0], label='3x', marker='o', color="tab:blue")
ax1.plot(all_data_fs[1], label='5x', marker='s', color="tab:orange")
ax1.plot(all_data_fs[2], label='10x', marker='D', color='tab:green')
ax1.plot(all_data_fs[3], label='30x', marker='d', color='tab:red')
ax1.set_xticks([1,5,10,15,20,25,30])
ax1.legend(title='Coverage')
ax1.set(xlabel="Derived allele frequency", ylabel="Number of sites")
ax1.text(-0.1, 1.05, 'A', transform=ax1.transAxes, size=12)


ax2.set_title("ANGSD Allele Frequency Spectra", fontsize=12)
ax2.plot(all_angsd_fs[0], label='3x', marker='o', color="tab:blue")
ax2.plot(all_angsd_fs[1], label='5x', marker='s', color="tab:orange")
ax2.plot(all_angsd_fs[2], label='10x', marker='D', color='tab:green')
ax2.plot(all_angsd_fs[3], label='30x', marker='d', color='tab:red')
ax2.legend(title='Coverage')
ax2.set_xticks([1,5,10,15,20,25,30,35,40])
ax2.set(xlabel="Derived allele frequency", ylabel="Number of sites")
ax2.text(-0.1, 1.05, 'B', transform=ax2.transAxes, size=12)

fig.tight_layout()
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)

plt.savefig(f"fs_gatk_angsd.png", dpi=150)
