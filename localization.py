import numpy as np
import seaborn as sns
import matplotlib.pylab as plt
# from scipy.misc import imread
from skimage.color import rgb2gray
import matplotlib.image as im
import os
import pickle
import pandas as pd
from scipy import stats


def relpath(filename):

    return os.path.join(os.path.dirname(__file__), filename)


def imheat(image, save_on=True, name1=None, name2=None, type=None):
    greyscale = rgb2gray(image / 255)
    if save_on:
        pickle.dump(greyscale,open("matfile.pkl", "wb+"))
    else:
        prev = pickle.load(open("matfile.pkl","rb"))
        temp = np.zeros(prev.shape)
        for i in range(temp.shape[0]):
            prev_indeces = np.where(prev[i] > 0)
            greyscale_indeces = np.where(greyscale[i] > 0)
            indeces = np.intersect1d(prev_indeces,greyscale_indeces)
            temp[i][indeces] = prev[i][indeces] + greyscale[i][indeces]
        greyscale = temp

    grid = np.zeros((16,16))
    for row in range(grid.shape[0]):
        for col in range(grid.shape[0]):
            cell_size = int(greyscale.shape[0]/grid.shape[0])
            start_row = row*cell_size
            start_col =col*cell_size
            grid[row, col] = np.sum(greyscale[start_row:start_row+cell_size, start_col:start_col+cell_size])
    grid = (grid - np.min(grid))/ (np.max(grid) - np.min(grid))
    ax = sns.heatmap(grid, linewidth=0.5)
    if not save_on:
        plt.title(name1 + "_" + name2 + "_" + type + "_coexpression_heatmap")
        plt.savefig("C:/Users/t-shbals/OneDrive - Microsoft/School/habiblab/images"+"/"+name1+"_"+name2+"_"+type+"_co_heat.png")
    else:
        plt.title(name1+"_" + type + "_heatmap")
        plt.savefig(
            "C:/Users/t-shbals/OneDrive - Microsoft/School/habiblab/images" + "/" + name1 + "_" + type + "_heat.png")
    plt.show()
    print("done" + str(save_on))


# image = im.imread(relpath("C3-MAX_5xfad new 3 slc38a1 gsn mfge8 dapi 40x-1.png"))
# imheat(image, True, "mfge8", "gsn", "AD")

def get_plaques_corners(x_centers, y_centers, r):
    rad = r/2
    top_left_x = x_centers - rad
    top_left_y = y_centers - rad
    return top_left_x.astype(np.int), top_left_y.astype(np.int)


def density(gfap,plaques, plaques_x, plaques_y, r):
    plaques_x_corners, plaques_y_corners = get_plaques_corners(plaques_x, plaques_y, r)
    sum = r*r
    with_p = []
    for i in range(len(plaques_x_corners)):
        outline_x = min(1023, plaques_x_corners[i]+r)
        outline_y = min(1023, plaques_y_corners[i]+r)
        temp = gfap[plaques_y_corners[i]:outline_y,plaques_x_corners[i]:outline_x]
        with_p.append(np.count_nonzero(temp)/sum)
    without_p = []
    for row in range(0,len(gfap), r):
        for col in range(0, len(gfap[0]), r):
            has_p = False
            for x in plaques_x:
                if col+r > x > col:
                    for y in plaques_y:
                        if row+r > y > row:
                            has_p = True
                            break
                if has_p:
                    break
            if not has_p:
                outline_x = min(1023, col + r)
                outline_y = min(1023, row + r)
                without_p.append(np.count_nonzero(gfap[row:outline_y,col:outline_x])/sum)
    print("with p:", with_p)
    print("without p:", without_p)
    with_p_avg = np.average(with_p)
    without_p_avg = np.average(without_p)
    print("with p avg:", with_p_avg)
    print("without p avg:", without_p_avg)
    stat, p = stats.ks_2samp(with_p, without_p)
    print("density statistic = ", stat)
    print("density p-value = ", p)
    # plt.bar(x = ["With plaque", "Without plaque"], height=[with_p_avg,without_p_avg])
    # sns.boxplot(y="Density",data=[with_p,without_p],order=["With plaque","Without plaque"])
    bo = plt.boxplot(x=[with_p,without_p],labels=["With plaque","Without plaque"],patch_artist=True)
    bo['boxes'][0].set(facecolor='crimson')
    bo['boxes'][1].set(facecolor='teal')
    plt.ylabel("Density")
    plt.title("Average Density of Gfap, p-value = " + str(round(p,7)))
    plt.savefig("C:/Users/t-shbals/PycharmProjects/habiblab/boxplotDensities.png")
    plt.show()


def get_pixels(im):
    pixels = []
    for row in range(len(im)):
        for col in range(len(im[0])):
            if im[row][col] > 0:
                pixels.append(np.array([row, col]))
    return np.array(pixels)


def find_min(curr_gfap_x, curr_gfap_y, plaques_pixels):
    dist = np.sqrt(np.power(plaques_pixels[:,0]-curr_gfap_x,2)+np.power(plaques_pixels[:,1]-curr_gfap_y,2))
    return np.min(dist)


def get_min_dists(gfap_pixels, plaques_pixels):
    gfap_min_dists = []
    for i in range(len(gfap_pixels)):
        gfap_min_dists.append(find_min(gfap_pixels[i][0], gfap_pixels[i][1], np.array(plaques_pixels)))
    return gfap_min_dists


def dist_vs_normal(gfap_min_dists):
    sns.distplot(gfap_min_dists, label='Minimum distance distribution', hist=False)
    normalDist = np.random.normal(loc=np.mean(gfap_min_dists), scale=50, size=len(gfap_min_dists))
    sns.distplot(normalDist, label="Normal distribution", hist=False)
    stat, p = stats.ks_2samp(gfap_min_dists, normalDist, mode="exact")
    print("normal stat = ", stat)
    print("normal p = ", p)
    plt.legend(fontsize="small")
    plt.xlabel("Distance")
    plt.ylabel("Percentage")
    plt.title("Distribution of minimum distances between GFAP and Plaques vs. Normal distribution\n"+"p-value = "+ str(round(p,325)), fontsize=10)
    plt.savefig("C:/Users/t-shbals/PycharmProjects/habiblab/distanceDistributionsVsNormalEmpty.png")
    plt.show()


def dist_vs_random(gfap_pixels, gfap_min_dists, num_of_plaques, radius, image_rows, image_cols):
    random_dists = []
    rand_obj_len = 2*radius
    for i in range(1):
        print(i)
        random_plaque_corners = np.random.randint(1023,size=(num_of_plaques,2))
        random_plaque_pixels = []
        for corner in random_plaque_corners:
            temp_row = corner[0]
            temp_col = corner[1]
            for r in range(temp_row, min(temp_row+rand_obj_len,image_rows)):
                for c in range(temp_col, min(temp_col+rand_obj_len,image_cols)):
                    random_plaque_pixels.append(np.array([r, c]))
        random_dists += get_min_dists(gfap_pixels,np.array(random_plaque_pixels))
    stat, p = stats.ks_2samp(gfap_min_dists, random_dists, mode="exact")
    print("random stat = ", stat)
    print("random p = ", p)
    sns.distplot(gfap_min_dists, label='Minimum distance distribution', hist=False)
    sns.distplot(random_dists, label="Random distribution", hist=False)
    plt.legend(fontsize="small")
    plt.xlabel("Distance")
    plt.ylabel("Percentage")
    plt.title("Distribution of minimum distances between GFAP and Plaques vs. Random distribution\n" +"p-value = "+ str(round(p,325)), fontsize=10)
    plt.savefig("C:/Users/t-shbals/PycharmProjects/habiblab/distanceDistributionsVsRandomEmpty.png")
    plt.show()



def calc_distances(gfap, plaques, num_of_plaques, r):
    gfap_pixels = get_pixels(gfap)
    plaques_pixels = get_pixels(plaques)
    # gfap_min_dists = []
    # for i in range(len(gfap_pixels)):
    #     gfap_min_dists.append(find_min(gfap_pixels[i][0],gfap_pixels[i][1], plaques_pixels))
    gfap_min_dists = get_min_dists(gfap_pixels, plaques_pixels)
    dist_vs_normal(gfap_min_dists)
    dist_vs_random(gfap_pixels, gfap_min_dists, num_of_plaques, r, len(plaques), len(plaques[0]))


def findR(plaques):
    r = []
    count = 0
    for row in range(len(plaques)):
        flag = False
        for col in range(len(plaques[0])):
            if plaques[row][col] > 0:
                flag = True
                count+=1
            else:
                if flag:
                    flag = False
                    if count > 5:
                        r.append(count)
                    count = 0
    print(np.max(np.array(r)))
    count = 0
    for col in range(len(plaques[0])):
        flag = False
        for row in range(len(plaques)):
            if plaques[row][col] > 0:
                flag = True
                count+=1
            else:
                if flag:
                    flag = False
                    # r = max(r, count)
                    if count > 9:
                        r.append(count)
                    count = 0
    print(np.max(np.array(r)))
    print(np.mean(np.array(r)))
    print(r)
    return int(round(np.average(np.array(r))))



def runGfapPlaqueAnalysis():
    gfap = im.imread(relpath("gfap_proteins.png")).astype(np.int)
    plaques = im.imread(relpath("plaques.png"))
    columns = pd.read_csv("MyExpt_IdentifyPrimaryObjectsRed.csv")
    columns = columns.fillna(0)
    relevantColumns = columns[["Location_Center_X", "Location_Center_Y"]]
    plaques_x = relevantColumns["Location_Center_X"].values
    plaques_y = relevantColumns["Location_Center_Y"].values
    num_of_plaques = len(plaques_x)
    r = findR(plaques)
    calc_distances(gfap, plaques, num_of_plaques, r)
    density(gfap,plaques,plaques_x,plaques_y,r)

runGfapPlaqueAnalysis()
