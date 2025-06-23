# Function to read Fortran unformatted binary files
# To use import the library using the following command:
# exec(open('/home/bilottos/Documents/PhD/Thesis_project/Programs/boutnp_sph_parallel/code/src/utils.py').read())

def readfbin(filePath):
    with open(filePath, 'rb') as f:
        _ = np.fromfile(f, dtype=np.int32, count=1)
        scalar = np.fromfile(f, dtype=np.int32, count=1)[0]
        _ = np.fromfile(f, dtype=np.int32, count=2)
        n1 = np.fromfile(f, dtype=np.int32, count=1)[0]
        n2 = np.fromfile(f, dtype=np.int32, count=1)[0]
        n3 = np.fromfile(f, dtype=np.int32, count=1)[0]
        _ = np.fromfile(f, dtype=np.int32, count=2)
        metadata = np.fromfile(f, dtype=np.float64, count=4)
        _ = np.fromfile(f, dtype=np.int32, count=2)
        modv = np.fromfile(f, dtype=np.float64, count=n1*n2*n3).reshape(n1,n2,n3, order='F')
        v_ph = np.fromfile(f, dtype=np.float64, count=n1*n2*n3).reshape(n1,n2,n3, order='F')
        v_rh = np.fromfile(f, dtype=np.float64, count=n1*n2*n3).reshape(n1,n2,n3, order='F')
        v_th = np.fromfile(f, dtype=np.float64, count=n1*n2*n3).reshape(n1,n2,n3, order='F')
    with open(filePath.parent / 'grid.xyz', 'rb') as f:
        _ = np.fromfile(f, dtype=np.int32, count=1)
        scalar = np.fromfile(f, dtype=np.int32, count=1)[0]
        _ = np.fromfile(f, dtype=np.int32, count=2)
        n1g = np.fromfile(f, dtype=np.int32, count=1)[0]
        n2g = np.fromfile(f, dtype=np.int32, count=1)[0]
        n3g = np.fromfile(f, dtype=np.int32, count=1)[0]
        _ = np.fromfile(f, dtype=np.int32, count=2)
        if (n1g==n1 and n2g==n2 and n3g==n3):
            xc = np.fromfile(f, dtype=np.float64, count=n1*n2*n3).reshape(n1,n2,n3, order='F')
            yc = np.fromfile(f, dtype=np.float64, count=n1*n2*n3).reshape(n1,n2,n3, order='F')
            zc = np.fromfile(f, dtype=np.float64, count=n1*n2*n3).reshape(n1,n2,n3, order='F')
        else:
            raise Exception("Field size incompatible with grid size")
    return {
        'n1':   n1,
        'n2':   n2,
        'n3':   n3,
        'x':    xc,
        'y':    yc,
        'z':    zc,
        'E':    modv,
        'vph':  v_ph,
        'vrh':  v_rh,
        'vth':  v_th
    }

def radial_rms(f, xc):
    n1 = f.shape[0]
    radius = xc[0,:,int((xc.shape[2]-1)/2)]
    sin_th = xc[0,0]/radius[0]
    sum_f = (f.sum(axis=0)*sin_th[np.newaxis,:]).sum(axis=1)
    avg_f = sum_f / (n1 * sin_th.sum())
    f2 = np.square(f - avg_f[np.newaxis,:,np.newaxis])
    sum_f2 = (f2.sum(axis=0)*sin_th[np.newaxis,:]).sum(axis=1)
    avg_f2 = sum_f2 / (n1 * sin_th.sum())
    return avg_f, np.sqrt(avg_f2)

def plotRadial(filePath, savePath=None):
    (n1,n2,n3,xc,yc,zc,E,vph,vrh,vth) = readfbin(filePath).values()
    radius = xc[0,:,int((xc.shape[2]-1)/2)]
    avg, rms = radial_rms(vrh, xc)
    plt.plot(radius, avg, c='xkcd:orange', ls='dashed', label='$\\bar{v}_r$')
    avg, rms = radial_rms(vth, xc)
    plt.plot(radius, avg, c='xkcd:sky blue', ls='dashed', label='$\\bar{v}_\\theta$')
    avg, rms = radial_rms(vph, xc)
    plt.plot(radius, avg, c='xkcd:pink', ls='dashed', label='$\\bar{v}_\\varphi$')
    avg, rms = radial_rms(vrh, xc)
    plt.plot(radius, rms, c='xkcd:orange', label='$\\sigma(v_r)$')
    avg, rms = radial_rms(vth, xc)
    plt.plot(radius, rms, c='xkcd:sky blue', label='$\\sigma(v_\\theta)$')
    avg, rms = radial_rms(vph, xc)
    plt.plot(radius, rms, c='xkcd:pink', label='$\\sigma(v_\\varphi)$')
    plt.legend(ncol=2)
    plt.xlabel('$r$')
    plt.ylabel('$v$')
    if savePath is None:
        plt.show()
    else:
        plt.title('Radial behaviour')
        plt.savefig(savePath)
        plt.close()
    
# plotRadial('spin_free/spbou20000.q')


def readspctr(filePath):
    with open(filePath,'rb') as f:
        _ = np.fromfile(f, dtype=np.int32, count=1)
        n3 = np.fromfile(f, dtype=np.int32, count=1)[0]
        n1mh = np.fromfile(f, dtype=np.int32, count=1)[0]
        _ = np.fromfile(f, dtype=np.int32, count=2)
        modes = np.fromfile(f, dtype=np.float64, count=n3*n1mh).reshape(n3,n1mh, order='F')
        _ = np.fromfile(f, dtype=np.int32, count=2)
        spctr = np.fromfile(f, dtype=np.float64, count=n3*n1mh).reshape(n3,n1mh, order='F')
        _ = np.fromfile(f, dtype=np.int32, count=2)
        time = np.fromfile(f, dtype=np.float64, count=1)[0]
        print('Imported arrays <modes> and <spctr> of size ' + str(n3) + 'x' + str(n1mh) + ', at time ' + str(time))
        return modes, spctr

def plotSpctr(spctr, modes, univScal=False, savePath=None):
    infernoCmap = plt.get_cmap('inferno')
    n3 = modes.shape[0]
    n1h = modes.shape[1]
    n3h = int(n3 / 2)
    if univScal:
        # Kolmogorov scaling
        eps = (modes**2 * spctr).sum(axis=1)
        eta = eps**(-1/4)
        for i in range(0, n3h):
            plt.loglog(modes[i,1:]*eta[i],(spctr[i,1:]/eps[i]**(1/4)+spctr[n3-1-i,1:]/eps[n3-i]**(1/4))/2, c=infernoCmap(i/n3h))
        plt.xlabel(r'$k\,\eta$')
        plt.ylabel(r'$E(k)\,\epsilon^{-1/4}$')
        plt.plot(modes[n3h,1:]*eta[n3h],1e-1*(modes[n3h,1:]*eta[n3h])**(-5/3))
    else:
        for i in range(0, n3h):
            plt.loglog(modes[i,1:],(spctr[i,1:]+spctr[n3-1-i,1:])/2, c=infernoCmap(i/n3h))
        plt.xlabel(r'$k$')
        plt.ylabel(r'$E(k)$')
        plt.plot(modes[n3h,1:],1e-1*modes[n3h,1:]**(-5/3))
    sm = plt.cm.ScalarMappable(cmap=infernoCmap, norm=plt.Normalize(vmin=0,vmax=np.pi))
    plt.colorbar(sm, label='polar angle [degrees]').ax.invert_yaxis()
    if savePath is None:
        plt.show()
    else:
        plt.title('Energy Spectrum')
        plt.savefig(savePath)

from scipy.ndimage import gaussian_filter1d
def fitMinMax(spctr, modes, i):
    xMin = []
    yMin = []
    xMax = []
    yMax = []
    plt.loglog(modes[i], spctr[i], label=f'$E(k)$ for $\\theta = {i/128*90}^\\circ$')
    plt.xlabel(r'$k$')
    plt.ylabel(r'$E(k)$')
    smooth_spctr = gaussian_filter1d(spctr[i,1:], 2)
    for j in range(smooth_spctr[1:].argmax(),smooth_spctr.size-2):
        if (smooth_spctr[j-1] > smooth_spctr[j] and smooth_spctr[j] < smooth_spctr[j+1]):
            xMin.append(modes[i,j])
            yMin.append(spctr[i,j])
        if (smooth_spctr[j-1] < smooth_spctr[j] and smooth_spctr[j] > smooth_spctr[j+1]):
            xMax.append(modes[i,j])
            yMax.append(spctr[i,j])
    plt.scatter(xMin, yMin, marker='s')
    mMin, cMin = np.polyfit(np.log(xMin), np.log(yMin), 1) # fit log(y) = m*log(x) + c
    plt.loglog(modes[i], np.exp(mMin*np.log(modes[i]) + cMin), linestyle='--', label=f'$m = {mMin}$')
    plt.scatter(xMax, yMax, marker='o')
    mMax, cMax = np.polyfit(np.log(xMax), np.log(yMax), 1) # fit log(y) = m*log(x) + c
    plt.loglog(modes[i], np.exp(mMax*np.log(modes[i]) + cMax), linestyle='--', label=f'$m = {mMax}$')
    plt.legend()
    plt.show()

# modes, spctr = readspctr(wdir/'noslip/ross/0.3/spctr080000.bin')
# plotSpctr(spctr, modes)
# fitMinMax(spctr, modes, 128)

from scipy.stats import binned_statistic
def avgSpctr(spctr, modes):
    n3 = modes.shape[0]
    n1q = int(modes.shape[1] / 2)
    binEdges = np.concatenate([np.linspace(1, n3/2, n1q)[:-1], np.logspace(np.log10(n3/2), np.log10(modes[2:-3].max()/2), n1q)])
    avgSpctr, binEdges, _ = binned_statistic(modes[2:-3,1:].flatten(), spctr[2:-3,1:].flatten(), statistic='mean', bins=binEdges)
    avgModes = (binEdges[1:] + binEdges[:-1])/2
    return avgModes, avgSpctr

# # Plot rescaled energy spectrum
# modes, spctr = readspctr(wdir/'free/23/spctr007500.bin'); plt.loglog(*avgSpctr(spctr*(9/2)**2, modes), label='2:3')
# modes, spctr = readspctr(wdir/'free/910/spctr046000.bin'); plt.loglog(*avgSpctr(spctr, modes), label='9:10')
# modes, spctr = readspctr(wdir/'free/1920/spctr070000.bin'); plt.loglog(*avgSpctr(spctr*(9/19)**2, modes), label='19:20')
# modes, spctr = readspctr(wdir/'free/4950/spctr066000.bin'); plt.loglog(*avgSpctr(spctr*(9/49)**2, modes), label='49:50')
# n3h = int(modes.shape[0] / 2)
# plt.plot(modes[n3h,1:],1e-1*modes[n3h,1:]**(-5/3))
# plt.xlabel(r'$k$')
# plt.ylabel(r'$E(k)$')
# plt.legend()
# plt.show()

# # Plot energy VS time
# data = np.genfromtxt(wdir / 'free/12/enerv.out').T;   plt.plot(data[0], data[1],       label='1:2')
# data = np.genfromtxt(wdir / 'free/34/enerv.out').T;   plt.plot(data[0], data[1]/3**2,  label='3:4')
# data = np.genfromtxt(wdir / 'free/89/enerv.out').T;   plt.plot(data[0], data[1]/8**2,  label='8:9')
# data = np.genfromtxt(wdir / 'free/1920/enerv.out').T; plt.plot(data[0], data[1]/19**2, label='19:20')
# data = np.genfromtxt(wdir / 'free/4950/enerv.out').T; plt.plot(data[0], data[1]/49**2, label='49:50')
# plt.xlabel('$t$')
# plt.ylabel('$E$')
# plt.legend()
# plt.show()