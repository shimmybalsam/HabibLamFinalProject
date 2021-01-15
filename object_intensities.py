import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr
# from lmfit.models import LorentzianModel
# from scipy.interpolate import make_interp_spline, BSpline

blue_intensities = []
red_intensities = []
green_intensities = []
cyon_intensities = []
all_intensities = []
cyon_counts = []
red_counts = []
green_counts = []
all_single_counts = []
all_single_counts_binary = []
all_cyon_red_counts = []
all_cyon_green_counts = []
all_green_red_counts = []
all_double_counts = []
all_double_freqs = []
all_triple_counts = []
# image_names = []
image_names = ['FAD','WT']
objects_count_per_image = []
colors = {0:'cyan', 1:'red', 2:'green'}
doubles = {0:'cyon_red',1:"cyon_green",2:"red_green"}
max_by_color = [0, 0, 0]

def clean_name(name):
    name = name.replace('.tif','')
    name = name.replace('.png','')
    return name


def get_foldchange(a,b):
    fold = a/b
    foldchange = round(math.log(fold,2),3)
    foldstr = "\n"+"Fold change of log(FAD/WT) = "+str(foldchange)+"\n"+"*log base 2"
    return foldstr

def find_threshold(cyon, red, green, image_label, i):
    # plt.hist(cyon, label='cyon', bins=20)
    # plt.hist(red, label='red', bins=20)
    # plt.hist(green, label='green', bins=20)
    # plt.legend()
    # plt.title(image_label)
    # plt.xlim(0,10)
    # plt.show()
    images = ['FAD1','FAD2','WT1','WT2']

    cyon.sort()
    red.sort()
    green.sort()
    plt.plot(cyon, label='cyon')
    plt.plot(red, label='red')
    plt.plot(green, label='green')
    plt.legend()
    plt.title(image_label)
    plt.ylim(top=10)
    plt.savefig('/Users/shimmy/PycharmProjects/per_object/findThreshold_'+images[i-1]+'.png')
    plt.show()

def count_setup():
    columns = pd.read_csv("MyExpt_IdentifySecondaryObjectsBlue.csv")
    columns = columns.fillna(0)
    relevantColumns = columns[["ImageNumber","ObjectNumber", "Metadata_File","Children_IdentifyPrimaryObjectsCyon_Count", "Children_IdentifyPrimaryObjectsGreen_Count", "Children_IdentifyPrimaryObjectsRed_Count", "Location_Center_X", "Location_Center_Y"]]
    numOfRows = relevantColumns.shape[0] - 1
    numOfImages = int(relevantColumns.iat[numOfRows,0])
    for i in range(1, numOfImages + 1):
        cyon_counts_binary = 0
        red_counts_binary = 0
        green_counts_binary = 0
        cyon_red_counts = 0
        cyon_green_counts = 0
        green_red_counts = 0
        triple_counts = 0
        relevantRows = relevantColumns.loc[relevantColumns['ImageNumber'] == i]
        image_label = relevantRows["Metadata_File"].values[0]
        # image_label = clean_name(image_label)
        # if image_label not in image_names:
        #     image_names.append(image_label)
        cyon_counts = relevantRows["Children_IdentifyPrimaryObjectsCyon_Count"].values
        red_counts = relevantRows["Children_IdentifyPrimaryObjectsRed_Count"].values
        green_counts = relevantRows["Children_IdentifyPrimaryObjectsGreen_Count"].values
        for j in range(len(cyon_counts)):
            if cyon_counts[j] > 0:
                cyon_counts_binary+=1
            if red_counts[j] > 0:
                red_counts_binary+=1
            if green_counts[j] > 0:
                green_counts_binary+=1
            if cyon_counts[j] > 0 and red_counts[j] > 0:
                cyon_red_counts+=1
            if cyon_counts[j] > 0 and green_counts[j] > 0:
                cyon_green_counts+=1
            if green_counts[j] > 0 and red_counts[j] > 0:
                green_red_counts+=1
            if green_counts[j] > 0 and red_counts[j] > 0 and cyon_counts[j] > 0:
                triple_counts+=1
        all_single_counts.append([cyon_counts,red_counts,green_counts])
        all_single_counts_binary.append([cyon_counts_binary,red_counts_binary,green_counts_binary])
        all_double_counts.append([cyon_red_counts,cyon_green_counts,green_red_counts])
        all_triple_counts.append(triple_counts)
        # find_threshold(cyon_counts,red_counts,green_counts,image_label, i)
        num_of_objects = int(relevantRows['ObjectNumber'].iloc[-1])
        objects_count_per_image.append(num_of_objects)
        max_by_color[0] = max(max(cyon_counts), max_by_color[0])
        max_by_color[1] = max(max(red_counts), max_by_color[1])
        max_by_color[2] = max(max(green_counts), max_by_color[2])
    # print(image_names)

def plot_objects_in_image():
    plt.bar(image_names,objects_count_per_image)
    plt.xticks(fontsize=6)
    plt.xlabel(get_foldchange(objects_count_per_image[0],objects_count_per_image[1]))
    plt.title("Number of cells/objects per image")
    plt.savefig('/Users/shimmy/PycharmProjects/per_object/numOfObjectsPerImage.png')
    plt.show()
    print(objects_count_per_image)

def check_count_averages():
    print("all double counts:")
    print(all_double_counts)
    print()
    print("objects_count_per_image:")
    print(objects_count_per_image)
    print("-----------------------------------------")
    print()
    triple_avg = []
    for color in range(3):
        single_avg = []
        single_binary_avg = []
        double_avg = []
        for i in range(len(all_single_counts)):
            single_avg.append(sum(all_single_counts[i][color])/objects_count_per_image[i])
            single_binary_avg.append(all_single_counts_binary[i][color]/objects_count_per_image[i])
            double_avg.append(all_double_counts[i][color]/objects_count_per_image[i])

        print(double_avg)
        plt.bar(x=['FAD','WT'],height = single_avg)
        # plt.bar(x=['FAD1','FAD2','WT1','WT2'],height = single_avg)
        plt.title('per_object_' + colors.get(color) + '_channel_average_single_count')
        plt.ylabel("Average")
        plt.xlabel(get_foldchange(single_avg[0],single_avg[1]))
        plt.savefig('/Users/shimmy/PycharmProjects/per_object/' + colors.get(color) + 'AverageSingleCount.png')
        plt.show()
        plt.bar(x=['FAD', 'WT'], height=single_binary_avg)
        # plt.bar(x=['FAD1', 'FAD2', 'WT1', 'WT2'], height=double_avg)
        plt.ylabel("Frequency")
        plt.xlabel(get_foldchange(single_binary_avg[0],single_binary_avg[1]))
        plt.title('per_object_' + colors.get(color) + '_channel_single_binary_count_frequency')
        plt.savefig('/Users/shimmy/PycharmProjects/per_object/' + colors.get(color) + 'SingleBinaryFreq.png')
        plt.show()
        plt.bar(x=['FAD', 'WT'], height=double_avg)
        # plt.bar(x=['FAD1', 'FAD2', 'WT1', 'WT2'], height=double_avg)
        plt.ylabel("Frequency")
        plt.xlabel(get_foldchange(double_avg[0],double_avg[1]))
        plt.title('per_object_' + doubles.get(color) + '_channels_double_count_frequency')
        plt.savefig('/Users/shimmy/PycharmProjects/per_object/' + doubles.get(color) + 'DoubleFreq.png')
        plt.show()
        # print(double_avg)
        # print(single_avg)
        # print("----------------------")
        # print()
    for i in range(len(all_triple_counts)):
        triple_avg.append(all_triple_counts[i]/objects_count_per_image[i])
    # print(all_triple_counts)
    # print(triple_avg)
    # print(objects_count_per_image)
    plt.bar(x=['FAD', 'WT'], height=triple_avg)
    # plt.bar(x=['FAD1', 'FAD2', 'WT1', 'WT2'], height=triple_avg)
    plt.title('per_object_' +'triple_channels_count_frequency')
    plt.ylabel("Frequency")
    plt.xlabel(get_foldchange(triple_avg[0],triple_avg[1]))
    plt.savefig('/Users/shimmy/PycharmProjects/per_object/tripleFreq.png')
    plt.show()


# def fit_func(x,y):
#     xnew = np.linspace(x.min(), x.max(), 300)
#     spl = make_interp_spline(x, y, k=1)
#     y_smooth = spl(xnew)
#     return xnew, y_smooth

def fit_func(x,a,b,c):
    return np.exp(-a*x+c)+b


def normalized_single_count():
    for color in range(3):
        normalized_across_images = {}
        normalized_arr = []
        x_axis = np.arange(max_by_color[color]+1)
        for i in range(len(all_single_counts)):
            vals, counts = np.unique(all_single_counts[i][color], return_counts=True)
            normalized = counts/np.sum(counts)
            normalized_padded = np.zeros(max_by_color[color]+1)
            normalized_padded[vals] = normalized
            normalized_arr.append(normalized_padded)
            normalized_across_images[image_names[i]] = normalized_padded

            # print("image"+str(i+1)+colors.get(color)+":")
            # print("vals:")
            # print(vals)
            # print("counts:")
            # print(counts)
            # print("normalized:")
            # print(normalized)
            # print("normalized_padded:")
            # print(normalized_padded)
            # print(np.sum(counts) == objects_count_per_image[i])
            # print("----------------------------------------")
            # print()
        df = pd.DataFrame(normalized_across_images, index=x_axis)
        df.plot.bar(rot=0)
        popt, pcov = curve_fit(fit_func, x_axis, normalized_arr[0])
        plt.plot(x_axis, fit_func(x_axis, *popt), label='FAD curve')
        popt, pcov = curve_fit(fit_func, x_axis, normalized_arr[1])
        plt.plot(x_axis, fit_func(x_axis, *popt), label='WT curve')
        plt.title('per_object_' + colors.get(color) + '_channel')
        plt.ylabel('Percent of cells')
        plt.xlabel('Cell count of channel')
        plt.legend()
        plt.savefig('/Users/shimmy/PycharmProjects/per_object/' + colors.get(color) + 'NormalizedBarsAcrossAllImages.png')
        plt.show()

        df = pd.DataFrame(normalized_across_images, index=x_axis)
        df.plot.bar(rot=0)
        # plt.legend(fontsize='small')
        popt, pcov = curve_fit(fit_func, x_axis, normalized_arr[0])
        plt.plot(x_axis, fit_func(x_axis, *popt), label='FAD curve')
        popt, pcov = curve_fit(fit_func, x_axis, normalized_arr[1])
        plt.plot(x_axis, fit_func(x_axis, *popt), label='WT curve')
        plt.title('per_object_' + colors.get(color) + '_channel')
        plt.xlim(left=6)
        plt.ylim(bottom= 0,top=0.02)
        plt.ylabel('Percent of cells')
        plt.xlabel('Cell count of channel')
        plt.legend()
        plt.savefig('/Users/shimmy/PycharmProjects/per_object/' + colors.get(color) + 'NormalizedBarsZoomAcrossAllImages.png')
        plt.show()

        df = pd.DataFrame(normalized_across_images, index=x_axis)
        df.plot(rot=0)
        popt, pcov = curve_fit(fit_func,x_axis,normalized_arr[0])
        plt.plot(x_axis,fit_func(x_axis,*popt),label='FAD curve')
        popt, pcov = curve_fit(fit_func, x_axis, normalized_arr[1])
        plt.plot(x_axis, fit_func(x_axis, *popt), label='WT curve')
        # df.plot(kind='area',rot=0)
        # plt.legend(fontsize='small')
        # x, y = fit_func(x_axis, normalized_arr[0])
        # plt.plot(x,y)
        # x, y = fit_func(x_axis, normalized_arr[1])
        # plt.plot(x, y)
        plt.title('per_object_' + colors.get(color) + '_channel')
        plt.xlim(left=5)
        plt.ylim(bottom = 0,top=0.02)
        plt.ylabel('Percent of cells')
        plt.xlabel('Cell count of channel')
        plt.legend()
        plt.savefig('/Users/shimmy/PycharmProjects/per_object/' + colors.get(color) + 'NormalizedPlotZoomAcrossAllImages.png')
        plt.show()


# def scatter_per_double():
#     markers = ["^", "*"]
#     m = 0
#     for im in all_single_counts:
#         x = np.arange(len(im[0]))
#         for i in range(3):
#             plt.scatter(x, im[i], c=colors.get(i), marker=markers[m])
#         plt.xlabel("Cell identity number")
#         plt.ylabel("Gene count")
#         plt.show()
#         m+=1
#
#     # plt.show()
#
# def scatter_double_diff():
#     markers = ["^", "*"]
#
#     m = 0
#     for im in all_single_counts:
#         x = np.arange(len(im[0]))
#         # plt.scatter(x, im[0], c=colors.get(0), marker=markers[m])
#         plt.scatter(x, np.array(im[0])-np.array(im[1]), c="cyan", marker=markers[m])
#         m += 1
#     plt.xlabel("Cell identity number")
#     plt.ylabel("Gene count")
#     plt.show()
#
#     m = 0
#     for im in all_single_counts:
#         x = np.arange(len(im[0]))
#         # plt.scatter(x, im[2], c=colors.get(2), marker=markers[m])
#         plt.scatter(x, np.array(im[2])-np.array(im[0]), c="g", marker=markers[m])
#         m += 1
#     plt.xlabel("Cell identity number")
#     plt.ylabel("Gene count")
#     plt.show()
#
#     m = 0
#     for im in all_single_counts:
#         x = np.arange(len(im[0]))
#         # plt.scatter(x, im[1], c=colors.get(1), marker=markers[m])
#         plt.scatter(x, np.array(im[1])-np.array(im[2]), c="r", marker=markers[m])
#         m += 1
#     plt.xlabel("Cell identity number")
#     plt.ylabel("Gene count")
#     plt.show()

def size_for_scatter(mat, tupple_dict):
    out_sizes = np.zeros(mat.shape[0])
    for i in range(len(mat)):
        tup = (mat[i][0], mat[i][1])
        out_sizes[i] = tupple_dict[tup]
    return out_sizes

def get_tupple_dict(unique, counts):
    d = {}
    for i in range(len(unique)):
        tup = (unique[i][0], unique[i][1])
        d[tup] = counts[i]
    return d

def scatter_per_double():
    markers = ["^", "*"]
    for i in range(3):
        genes = {0:"mfge8",1:"myoc",2:"slc38a1"}
        # genes = {0:"mfge8",1:"myoc",2:"osmr"}
        # genes = {0:"vim",1:"gfap",2:"slc38a1"}
        # genes = {0:"mfge8",1:"gsn",2:"slc38a1"}
        for m, im in enumerate(all_single_counts):
            # a = np.array(im[i])
            # b = np.array(im[(i + 1) % 3])
            # ind1 = np.where(a > 0)
            # ind2 = np.where(b >0)
            # indeces = np.intersect1d(ind1,ind2)
            # x = a[indeces]
            # y = b[indeces]
            x = np.array(im[i])
            y = np.array(im[(i + 1) % 3])
            temp = np.array([x,y])
            unique, counts = np.unique(temp,axis=1,return_counts=True)
            tupple_dict = get_tupple_dict(unique.T,counts)
            z = size_for_scatter(temp.T,tupple_dict).astype(np.int32)
            # print(z)
            model = LinearRegression().fit(x.reshape(-1, 1), y.reshape(-1, 1), sample_weight=z)
            r2 = round(model.score(x.reshape(-1, 1),y.reshape(-1, 1),sample_weight=z),2)
            r = round(pearsonr(x,y)[0],2)
            # r_check = np.corrcoef(x,y)
            # print(r)
            # print(r_check)

            # y_pred = model.intercept_ + model.coef_ * x.reshape(-1, 1)
            y_pred = model.predict(x.reshape(-1,1))

            gene1 = genes.get(i)
            gene2 = genes.get((i + 1) % 3)
            plt.scatter(x, y, c=colors.get(i), s=z, marker=markers[m])
            plt.plot(x, y_pred, linewidth=1)
            plt.title(image_names[m]+" correlation of "+gene1+" and "+gene2 + r", $r$ = "+str(r))
            plt.xlabel(gene1)
            plt.ylabel(gene2)
            plt.savefig("C:/Users/t-shbals/OneDrive - Microsoft/school/habiblab/images/"+image_names[m]+"_corr_"+gene1+"_"+gene2)
            plt.show()

    # plt.show()

def scatter_double_diff():
    markers = ["^", "*"]

    m = 0
    for im in all_single_counts:
        x = np.arange(len(im[0]))
        # plt.scatter(x, im[0], c=colors.get(0), marker=markers[m])
        plt.scatter(np.array(im[1]), np.array(im[0])-np.array(im[1]), c="cyan", marker=markers[m])
        m += 1
    plt.xlabel("Red")
    plt.ylabel("Cyan - Red")
    plt.show()

    m = 0
    for im in all_single_counts:
        x = np.arange(len(im[0]))
        # plt.scatter(x, im[2], c=colors.get(2), marker=markers[m])
        plt.scatter(np.array(im[0]), np.array(im[2])-np.array(im[0]), c="g", marker=markers[m])
        m += 1
    plt.xlabel("Cyan")
    plt.ylabel("Green - Cyan")
    plt.show()

    m = 0
    for im in all_single_counts:
        x = np.arange(len(im[0]))
        # plt.scatter(x, im[1], c=colors.get(1), marker=markers[m])
        plt.scatter(np.array(im[2]), np.array(im[1])-np.array(im[2]), c="r", marker=markers[m])
        m += 1
    plt.xlabel("Green")
    plt.ylabel("Red  - Green")
    plt.show()


def count_hist_per_channel():
    for color in range(3):
        temp = []
        for i in range(len(all_single_counts)):

        #better way: make array of histogram which counts how many cells have each count and then divide that amount by number of cells:
        # as in: num of cells with value X / num of cells in image
        # try to use numpy functions, such as np.unique and make param return_counts=true.

            #side by side:
            temp.append(all_single_counts[i][color])
            # temp.append([j / objects_count_per_image[i] for j in all_single_counts[i][color]])

            #density:
            sns.distplot(temp[i], label=image_names[i], hist=False, kde=True)
        # print(temp[i])
        # print(all_single_counts[i][color])
        #density:

        # plt.legend(fontsize='small')
        # plt.title('per_object_' + colors.get(color) + '_channel')
        # plt.ylabel('Percent of Cells')
        # plt.xlabel('Cell count')
        # plt.savefig('/Users/shimmybalsam/PycharmProjects/per_object/' + colors.get(color) + 'DensityAcrossAllImages.png')
        # plt.show()

        # plt.xlim(-2,10)
        plt.legend(fontsize='small')
        plt.title('per_object_' + colors.get(color) + '_channel')
        plt.ylabel('Percent of Cells')
        plt.xlabel('Cell count')

        # plt.subplot(212)
        plt.xlim(left=10)
        plt.ylim(top=0.01)
        # plt.legend(fontsize='small')
        # # plt.title('per_object_' + colors.get(color) + '_channel')
        # plt.ylabel('Percent of Cells')
        # plt.xlabel('Cell count')
        plt.savefig('/Users/shimmy/PycharmProjects/per_object/' + colors.get(color) + 'DensityZoom2AcrossAllImages.png')
        plt.show()

        #side by side:
        plt.figure()
        plt.subplot(211)
        plt.hist(temp, label=image_names, bins=15)
        plt.legend(fontsize='small')
        plt.title('per_object_'+ colors.get(color)+'_channel')
        plt.ylabel('Amount of cells')
        plt.xlabel('Cell count of channel')

        plt.subplot(212)
        plt.hist(temp, label=image_names, bins=15)
        plt.legend(fontsize='small')
        # plt.title('per_object_' + colors.get(color) + '_channel')
        # plt.xlim(left=2)
        plt.xlim(left=0.01)
        plt.ylim(top=25)
        plt.ylabel('Amount of cells')
        plt.xlabel('Cell count of channel')
        plt.savefig('/Users/shimmy/PycharmProjects/per_object/' + colors.get(color) + 'HistAcrossAllImages.png')
        plt.show()
    # print(temp)


def count_sorted_per_channel():
    for color in range(3):
        for i in range(len(all_single_counts)):
            # all_single_counts[i][color].sort()
            # plt.plot(all_single_counts[i][color], label=image_names[i])
            normalized = [j / objects_count_per_image[i] for j in all_single_counts[i][color]]
            normalized.sort()
            # print(normalized)
            plt.plot(normalized, label=image_names[i])

        plt.legend(fontsize='small')
        plt.title('per_object_'+ colors.get(color)+'_channel')
        plt.xlabel('Cell Index')
        plt.ylabel('Channel count per cell')
        # plt.yscale('log')
        # plt.ylim(bottom=0)
        plt.savefig('/Users/shimmy/PycharmProjects/per_object/' + colors.get(color) + 'SortedPlotAcrossAllImages.png')
        plt.show()


# def count_hist_per_channel():
#
#
# def count_density_per_channel():
#
# def count_doubles_coexpression():
#
# def count_triple_coexpression():


def intensity_setup():
    for color in ['Red', 'Green', 'Cyon']:
        primaryColumns = pd.read_csv("MyExpt_IdentifyPrimaryObject" + color + ".csv")
        primaryColumns = primaryColumns.fillna(0)
        primaryNumOfRows = primaryColumns.shape[0] - 1

        relateColumns = pd.read_csv("MyExpt_relateObjectsBlue" + color + ".csv")
        relateColumns = primaryColumns.fillna(0)
        relateNumOfRows = primaryColumns.shape[0] - 1


def main():
    count_setup()
    scatter_per_double()
    # scatter_double_diff()
    # count_hist_per_channel()
    # count_sorted_per_channel()
    # normalized_single_count()
    # check_count_averages()
    # plot_objects_in_image()
main()


