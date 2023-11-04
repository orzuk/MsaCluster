from protein_utils import *
import nglview as nv
import py3Dmol
from pyvirtualdisplay import Display
from IPython.display import display, Image
import requests
from PIL import Image
from io import BytesIO
from MSA_Clust import match_predicted_and_true_contact_maps


# Plot structure using nglview
def plot_pdb_struct(pdb_file, output_3d_image_file = []):
#    view = nv.show_file(pdb_file)
#    view

    with open(pdb_file, "r") as f:
        pdb_data = f.read()

    # Create a 3Dmol view
    viewer = py3Dmol.view(width=400, height=400)

    # Add the protein structure data to the viewer
    viewer.addModel(pdb_data, "pdb")

    # Style the visualization (optional)
    viewer.setStyle({"cartoon": {"color": "spectrum"}})

    # Zoom to fit the structure
    viewer.zoomTo()

    output_file = f"protein_structure.png"
    with open(output_file, "wb") as f:
        f.write(viewer.toImage("png"))
    print(f"Saved {output_file}")

#    display(viewer)
#    output_file = "protein_structure.png"  # Replace with your desired output file name

#    with Display():
#        viewer.savefig(output_file)

    # Show the 3D visualization
#    viewer.show()

    if len(output_3d_image_file) == 0:
        output_3d_image_file = "Pipeline/Results/Figures/3d_struct/" + os.path.basename(pdb_file) + ".png"


#    screenshot = viewer.png()
#    with open(output_3d_image_file, 'wb') as f:
#        f.write(screenshot)

#    viewer.png(output_3d_image_file) # , width=400, height=400)

    return 0


# Plot multiple contacts and predictions together
def plot_array_contacts_and_predictions(predictions, contacts, save_file=[]):
    n_pred = len(predictions)
    n_row = math.ceil(math.sqrt(n_pred))  # *2
    if n_row * (n_row - 1) >= n_pred:  # *2
        n_col = n_row - 1
    else:
        n_col = n_row
    PDB_IDS = predictions.keys()  # [p[name] for p in predictions]
    contacts_ids = contacts.keys()
    fig, axes = plt.subplots(figsize=(18, 18), nrows=n_row, ncols=n_col, layout="compressed")
    #    print("Num cmaps: " + str(n_pred))
    #    print(axes.shape)
    #    fig, axes = plt.subplots(figsize=(18, 6), ncols=n_pred)
    ctr = 0
    #    for ax, name in zip(axes, PDB_IDS):
#    print("n_col=" + str(n_col))
#    print("n_row=" + str(n_row))
#    print(PDB_IDS)
#    print("Contact lens:" + str(len(contacts)))
    for name in PDB_IDS:  # loop over predictions
        if n_col == 1:
            ax = axes[ctr]
        else:
            ax = axes[ctr // n_col, ctr % n_col]
        ctr = ctr + 1
#        print("Plotting prediction: " + name)  # + " -> true: " + true_name)
        plot_foldswitch_contacts_and_predictions(
            predictions[name], contacts, ax=ax, title=name, show_legend = ctr == 1)

    ##        for true_name in contacts_ids: # loop over two folds
    ##            print("Plotting prediction: " + name + " -> true: " + true_name)
    ##           plot_contacts_and_predictions(
    ##                predictions[name], contacts[true_name], ax=ax, title = name)
    #            prediction, target, ax=ax, title = lambda prec: f"{name}: Long Range P@L: {100 * prec:0.1f}")
    if len(save_file) > 0:  # save and close plot (enable automatic saving of multiple plots)
        plt.savefig(save_file + '.png')
        print("Save cmap fig: " + save_file + '.png')
    else:
        plt.show()


"""Adapted from: https://github.com/rmrao/evo/blob/main/evo/visualize.py"""


def plot_contacts_and_predictions(
        predictions: Union[torch.Tensor, np.ndarray],
        contacts: Union[torch.Tensor, np.ndarray],
        ax: Optional[mpl.axes.Axes] = None,
        # artists: Optional[ContactAndPredictionArtists] = None,
        cmap: str = "Blues",
        ms: float = 1,
        title: Union[bool, str, Callable[[float], str]] = True,
        animated: bool = False,
) -> None:
    if isinstance(predictions, torch.Tensor):
        predictions = predictions.detach().cpu().numpy()
    if isinstance(contacts, torch.Tensor):
        contacts = contacts.detach().cpu().numpy()
    if ax is None:
        ax = plt.gca()

    seqlen = contacts.shape[0]
    relative_distance = np.add.outer(-np.arange(seqlen), np.arange(seqlen))
    bottom_mask = relative_distance < 0
    masked_image = np.ma.masked_where(bottom_mask, predictions)
    invalid_mask = np.abs(np.add.outer(np.arange(seqlen), -np.arange(seqlen))) < 6
    predictions = predictions.copy()
    predictions[invalid_mask] = float("-inf")

    topl_val = np.sort(predictions.reshape(-1))[-seqlen]
    pred_contacts = predictions >= topl_val
    true_positives = contacts & pred_contacts & ~bottom_mask
    false_positives = ~contacts & pred_contacts & ~bottom_mask
    other_contacts = contacts & ~pred_contacts & ~bottom_mask

    if isinstance(title, str):
        title_text: Optional[str] = title
    elif title:
        long_range_pl = compute_precisions(predictions, contacts, minsep=24)[
            "P@L"
        ].item()
        if callable(title):
            title_text = title(long_range_pl)
        else:
            title_text = f"Long Range P@L: {100 * long_range_pl:0.1f}"
    else:
        title_text = None

    img = ax.imshow(masked_image, cmap=cmap, animated=animated)  # Show main image
    oc = ax.plot(*np.where(other_contacts), "o", c="grey", ms=ms, label="other")[0]
    fp = ax.plot(*np.where(false_positives), "o", c="r", ms=ms, label="FP")[0]
    tp = ax.plot(*np.where(true_positives), "o", c="b", ms=ms, label="TP")[0]
    ti = ax.set_title(title_text) if title_text is not None else None
    # artists = ContactAndPredictionArtists(img, oc, fp, tp, ti)

    # Show second structure here!

    ax.legend(loc="upper left")
    ax.axis("square")
    ax.set_xlim([0, seqlen])
    ax.set_ylim([0, seqlen])
    save_flag = False  # add as input
    if save_flag:
        plt.savefig('%s.pdf' % title, bbox_inches='tight')


# Plot contacts and predictions for TWO folds !!!
def plot_foldswitch_contacts_and_predictions(
        predictions: Union[torch.Tensor, np.ndarray],
        contacts: Union[torch.Tensor, np.ndarray],
        ax: Optional[mpl.axes.Axes] = None,
        # artists: Optional[ContactAndPredictionArtists] = None,
        cmap: str = "gray_r",  # "Blues",
        ms: float = 5,
        title: Union[bool, str, Callable[[float], str]] = True,
        animated: bool = False,
        show_legend: bool = False,
) -> None:
    fold_ids = list(contacts.keys())
    if isinstance(predictions, torch.Tensor):
        predictions = predictions.detach().cpu().numpy()
    for fold in fold_ids:
        if isinstance(contacts[fold], torch.Tensor):
            contacts[fold] = contacts[fold].detach().cpu().numpy()
    if ax is None:
        ax = plt.gca()

    if len(fold_ids) == 1: # same PDB ID, duplicate
        fold_ids[1] = fold_ids[0]
    seqlen = contacts[fold].shape[0]
#    print(seqlen)
#    for fold in fold_ids:
#        print(contacts[fold].shape)
#    print(fold_ids)
    relative_distance = np.add.outer(-np.arange(seqlen), np.arange(seqlen))
    top_bottom_mask = {fold_ids[0]: relative_distance < 0, fold_ids[1]: relative_distance > 0}
    #    masked_image = np.ma.masked_where(bottom_mask, predictions)
    masked_image = np.ma.masked_where(top_bottom_mask[list(fold_ids)[0]], predictions)
    invalid_mask = np.abs(np.add.outer(np.arange(seqlen), -np.arange(seqlen))) < 6
    predictions = predictions.copy()
    predictions[invalid_mask] = float("-inf")

##    contacts_united = (contacts[fold_ids[0]] + contacts[fold_ids[1]])  # 0: no contact, 1: contact in one, 2: contact in both
##    for fold in fold_ids:
##        contacts_united[np.where(contacts[fold] & (contacts_united == 1) & top_bottom_mask[fold])] = 0
    # Flip contact in one and both:
 ##   cc = copy.deepcopy(contacts_united)
##    contacts_united[cc == 1] = 2
##    contacts_united[cc == 2] = 1

    contacts_united = match_predicted_and_true_contact_maps(predictions, contacts)

    topl_val = np.sort(predictions.reshape(-1))[-seqlen]
    pred_contacts = predictions >= topl_val
    true_positives, false_positives, other_contacts = {}, {}, {}  # [None]*2, [None]*2, [None]*2

    for fold in fold_ids:
        #        print(fold)
        #        print(true_positives[fold])
        #        print(contacts[fold])
        #        print("Top-Bottom")
        #       print(top_bottom_mask[fold])
        true_positives[fold] = contacts[fold] & pred_contacts & top_bottom_mask[fold]
        false_positives[fold] = ~contacts[fold] & pred_contacts & top_bottom_mask[fold]
        other_contacts[fold] = contacts[fold] & ~pred_contacts & top_bottom_mask[fold]

    if isinstance(title, str):
        title_text: Optional[str] = title
    elif title:
        long_range_pl = compute_precisions(predictions, contacts, minsep=24)[
            "P@L"
        ].item()
        if callable(title):
            title_text = title(long_range_pl)
        else:
            title_text = f"Long Range P@L: {100 * long_range_pl:0.1f}"
    else:
        title_text = None

    #    img = ax.imshow(masked_image, cmap=cmap, animated=animated)  # Show main image
    img = ax.imshow(contacts_united, cmap=cmap, animated=animated)  # Show main image
#    for fold in fold_ids:
#        oc = ax.plot(*np.where(other_contacts[fold]), "o", c="grey", ms=ms, label="other")[0]
    ms = ms * 50 / seqlen
#    print("ms: " + str(ms))
    fp = ax.plot(*np.where(false_positives[fold_ids[0]]), "o", c="r", ms=ms, label="FP")[0]
    tp = ax.plot(*np.where(true_positives[fold_ids[0]]), "o", c="b", ms=ms, label="TP")[0]
    fp = ax.plot(*np.where(false_positives[fold_ids[1]]), "o", c="r", ms=ms)[0]
    tp = ax.plot(*np.where(true_positives[fold_ids[1]]), "o", c="b", ms=ms)[0]
    ti = ax.set_title(title_text) if title_text is not None else None
    # artists = ContactAndPredictionArtists(img, oc, fp, tp, ti)

    # Show second structure here!
    if show_legend:
        ax.legend(loc="upper left")
    ax.axis("square")
    ax.set_xlabel(fold_ids[0], fontsize=14)
    ax.set_ylabel(fold_ids[1], fontsize=14)
    ax.set_xlim([0, seqlen])
    ax.set_ylim([0, seqlen])
    save_flag = False  # add as input
    if save_flag:
        plt.savefig('%s.pdf' % title, bbox_inches='tight')