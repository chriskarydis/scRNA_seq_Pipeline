# main.py
from io import BytesIO
import streamlit as st
import scanpy as sc
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse
import scanpy.external as sce
import plotly.graph_objects as go
import seaborn as sns
import warnings
import tempfile
import os
warnings.filterwarnings("ignore", message="set_ticklabels()", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


# Ρύθμιση σελίδας
st.set_page_config(page_title="scRNA-seq Pipeline", layout="wide")

# Τίτλος
st.title("🔬 Ανάλυση scRNA-seq δεδομένων")

# Tabs
tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "📁 Δεδομένα",
    "⚙️ Προεπεξεργασία",
    "📊 Ανάλυση",
    "🧬 Γονιδιακή Ανάλυση",
    "📈 Εξαγωγή Αποτελεσμάτων",
    "👥 Ομάδα"
])

# Tab 1: Προεπισκόπηση
with tab1:
    st.header("📁 Φόρτωση & Προεπισκόπηση Δεδομένων")

    uploaded_file = st.file_uploader("📤 Ανέβασε αρχείο τύπου `.h5ad`", type=["h5ad"])

    if uploaded_file is not None:
        try:
            adata = sc.read_h5ad(uploaded_file)
            st.session_state["adata"] = adata

            # Αν δεν υπάρχει ήδη, αρχικοποιούμε αντίγραφα
            if "adata_pre" not in st.session_state:
                st.session_state["adata_pre"] = adata.copy()
                st.session_state["preprocessing_done"] = False

            st.success("✅ Το αρχείο φορτώθηκε επιτυχώς!")

            # Γρήγορη σύνοψη
            col1, col2 = st.columns(2)
            col1.metric("🔬 Αριθμός Κυττάρων", f"{adata.n_obs:,}")
            col2.metric("🧬 Αριθμός Γονιδίων", f"{adata.n_vars:,}")

            # Expanders για obs και var
            with st.expander("🔍 Προεπισκόπηση `obs` (μεταδεδομένα κυττάρων)"):
                if adata.obs.shape[1] > 0:
                    obs_page = st.number_input("📄 Σελίδα metadata (10/σελ)", min_value=1, max_value=(len(adata.obs) - 1) // 10 + 1, step=1)
                    st.dataframe(adata.obs.iloc[(obs_page - 1) * 10 : obs_page * 10], use_container_width=True)
                else:
                    st.info("ℹ️ Δεν υπάρχουν μεταδεδομένα για τα κύτταρα.")

            with st.expander("🧬 Προεπισκόπηση `var` (γονίδια)"):
                if adata.var_names is not None and len(adata.var_names) > 0:
                    var_page = st.number_input("📄 Σελίδα γονιδίων (10/σελ)", min_value=1, max_value=(len(adata.var_names) - 1) // 10 + 1, step=1)
                    genes_df = pd.DataFrame(adata.var_names, columns=["Gene Names"])
                    st.dataframe(genes_df.iloc[(var_page - 1) * 10 : var_page * 10], use_container_width=True)
                else:
                    st.info("ℹ️ Δεν υπάρχουν ονόματα γονιδίων.")


        except Exception as e:
            st.error(f"❌ Σφάλμα κατά τη φόρτωση του αρχείου: {str(e)}")
    else:
        st.info("🔁 Ανέβασε αρχείο για να ξεκινήσεις.")

    st.markdown("---")
    st.info("➡️ Κάνε scroll επάνω και επίλεξε το επόμενο tab για να συνεχίσεις με την προεπεξεργασία.")


# Tab 2: Προεπεξεργασία
with tab2:
    st.header("⚙️ Προεπεξεργασία Δεδομένων")
    if "adata" not in st.session_state:
        st.warning("❗ Αρχικά ανέβασε αρχείο `.h5ad` στο πρώτο tab.")
    else:
        ad = st.session_state["adata"].copy()

        st.subheader("🔧 Επιλογές Προεπεξεργασίας")

        col1, col2 = st.columns(2)
        with col1:
            min_genes = st.slider("🧬 Ελάχιστα γονίδια ανά κύτταρο", min_value=0, max_value=1000, value=600, step=10)
            do_filter_cells = st.checkbox("🧪 Φιλτράρισμα Κυττάρων")

        with col2:
            min_cells = st.slider("🔬 Ελάχιστα κύτταρα ανά γονίδιο", min_value=0, max_value=50, value=3, step=1)
            do_filter_genes = st.checkbox("🧬 Φιλτράρισμα Γονιδίων")

        do_remove_mt = st.checkbox("❌ Αφαίρεση MT-/ERCC γονιδίων", value=True)
        do_normalize_log1p = st.checkbox("⚖️ Κανονικοποίηση + log1p", value=True)
        do_hvg = st.checkbox("🔬 Επιλογή HVGs", value=True)
        do_scale = st.checkbox("📏 Scaling", value=True)

        if st.button("🚀 Εκτέλεση Προεπεξεργασίας"):
            with st.status("🔄 Εκτελείται Προεπεξεργασία...", expanded=True) as status:
                try:
                    if do_filter_cells:
                        with st.spinner("🔬 Φιλτράρισμα Κυττάρων..."):
                            before = ad.shape[0]
                            sc.pp.filter_cells(ad, min_genes=min_genes)
                            after = ad.shape[0]
                            st.success(f"✅ Αφαιρέθηκαν {before - after} κύτταρα (έμειναν {after})")


                    if do_filter_genes:
                        with st.spinner("🧬 Φιλτράρισμα Γονιδίων..."):
                            before = ad.shape[1]
                            sc.pp.filter_genes(ad, min_cells=min_cells)
                            after = ad.shape[1]
                            st.success(f"✅ Αφαιρέθηκαν {before - after} γονίδια (έμειναν {after})")


                    if do_remove_mt:
                        with st.spinner("🧫 Αφαίρεση MT-/ERCC γονιδίων..."):
                            before = ad.shape[1]
                            genes_to_remove = [g for g in ad.var_names if g.upper().startswith(("MT-", "ERCC"))]
                            ad = ad[:, [g for g in ad.var_names if g not in genes_to_remove]]
                            after = ad.shape[1]
                            st.success(f"✅ Αφαιρέθηκαν {before - after} MT-/ERCC γονίδια (έμειναν {after})")


                    if do_normalize_log1p:
                        with st.spinner("⚖️ Κανονικοποίηση και log1p..."):
                            sc.pp.normalize_total(ad, target_sum=1e4)
                            sc.pp.log1p(ad)
                            st.success("✅ Κανονικοποίηση και log1p εφαρμόστηκαν")

                    if do_hvg:
                        with st.spinner("🔍 Επιλογή HVGs..."):
                            if scipy.sparse.issparse(ad.X):
                                ad.X = ad.X.toarray()
                            ad.X = np.nan_to_num(ad.X, nan=0.0, posinf=0.0, neginf=0.0)
                            sc.pp.highly_variable_genes(ad, min_mean=0.0125, max_mean=3, min_disp=0.5)
                            ad.raw = ad
                            before = ad.shape[1]
                            ad = ad[:, ad.var.highly_variable]
                            after = ad.shape[1]
                            st.success(f"✅ Επιλέχθηκαν {after} HVGs (από {before})")


                    if do_scale:
                        with st.spinner("📏 Εφαρμογή Scaling..."):
                            sc.pp.scale(ad, max_value=10)
                            st.success("✅ Scaling εφαρμόστηκε")

                    st.session_state["adata_pre"] = ad.copy()
                    st.session_state["preprocessing_done"] = True
                    status.update(label="🎉 Η προεπεξεργασία ολοκληρώθηκε", state="complete")

                except Exception as e:
                    st.error(f"❌ Σφάλμα: {str(e)}")
                    status.update(label="❌ Σφάλμα κατά την προεπεξεργασία", state="error")

        if st.button("🔁 Επαναφορά Αρχικών Δεδομένων"):
            with st.status("🔄 Επαναφορά αρχικών δεδομένων...", expanded=True) as status:
                with st.spinner("🔁 Γίνεται επαναφορά..."):
                    st.session_state["adata_pre"] = st.session_state["adata"].copy()
                    st.session_state["preprocessing_done"] = False
                    st.success("✅ Τα δεδομένα επανήλθαν στην αρχική τους μορφή")
                status.update(label="✅ Επαναφορά ολοκληρώθηκε", state="complete")

        # Υπολογισμός UMAP (πριν την ανάλυση)
        with st.spinner("📍 Υπολογισμός UMAP για οπτικοποίηση πριν την ανάλυση..."):
            try:
                sc.pp.pca(ad, n_comps=30)
                sc.pp.neighbors(ad, n_neighbors=10)
                sc.tl.umap(ad, n_components=3)
                st.session_state["adata_pre_umap"] = ad.copy()
                st.success("✅ Υπολογίστηκε UMAP (πριν την ανάλυση). Παρουσίαση στο επόμενο tab.")
            except Exception as e:
                st.warning(f"⚠️ Σφάλμα στον υπολογισμό του αρχικού UMAP: {str(e)}")



    st.markdown("---")
    st.info("➡️ Κάνε scroll επάνω και επίλεξε το επόμενο tab για να συνεχίσεις με την ανάλυση.")


# Tab 3: PCA, Clustering, UMAP & Harmony
with tab3:
    st.header("📊 Ανάλυση: PCA, Clustering, UMAP, Harmony")
    if "adata" not in st.session_state:
        st.warning("❗ Αρχικά ανέβασε αρχείο `.h5ad` στο πρώτο tab.")
    elif "adata_pre" not in st.session_state:
            st.warning("❗ Πρώτα ολοκλήρωσε την προεπεξεργασία.")
    else:

        adata = st.session_state["adata_pre"].copy()

        # Επιλογή: Χρήση Harmony;
        use_harmony = st.checkbox("🔧 Χρήση Harmony για διόρθωση batch effect", value=False)

        # Παράμετροι
        n_comps = st.slider("🔢 Αριθμός PCA components", min_value=2, max_value=100, value=30)
        n_neighbors = st.slider("👥 Γείτονες για Clustering", min_value=2, max_value=50, value=10)
        resolution = st.slider("📏 Leiden resolution", min_value=0.1, max_value=2.0, value=0.5, step=0.1)

        projection = st.radio("📐 Προβολή UMAP", options=["2D", "3D"], index=0, horizontal=True, key="umap_proj")
        n_umap_components = 3 if projection == "3D" else 2


        if st.button("🚀 Εκτέλεση Ανάλυσης"):
            with st.status("🔄 Εκτελείται Ανάλυση...", expanded=True) as status:
                try:
                    with st.spinner("🔢 Υπολογισμός PCA..."):
                        sc.pp.pca(adata, n_comps=n_comps)
                        st.success(f"✅ PCA ολοκληρώθηκε με {n_comps} components")
                        st.info("ℹ️ Το PCA μειώνει τη διάσταση των δεδομένων, διατηρώντας τις πιο σημαντικές συνιστώσες.")

                    if use_harmony:
                        with st.spinner("🔧 Harmony integration..."):
                            sce.pp.harmony_integrate(adata, key='batch')
                            st.success("✅ Εφαρμόστηκε Harmony για διόρθωση batch effect")
                            st.info("ℹ️ Το Harmony μειώνει τη μεταβλητότητα μεταξύ batch, ενοποιώντας τα δεδομένα.")

                            sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=n_neighbors)
                            st.success(f"✅ Υπολογίστηκαν γείτονες ({n_neighbors}) με βάση το PCA-Harmony")

                            sc.tl.leiden(adata, resolution=resolution, flavor="igraph", directed=False, random_state=0)
                            st.success(f"✅ Εκτελέστηκε Clustering με Leiden (resolution={resolution})")
                            st.info("ℹ️ Το Leiden εντοπίζει ομάδες κυττάρων βάσει γειτνίασης.")

                            sc.tl.umap(adata, n_components=n_umap_components)
                            st.success("✅ UMAP ολοκληρώθηκε πάνω στο PCA-Harmony")
                            st.info("ℹ️ Το UMAP δημιουργεί δισδιάστατη προβολή για εύκολη οπτικοποίηση.")
                    else:
                        with st.spinner("👥 Clustering & UMAP χωρίς Harmony..."):
                            sc.pp.neighbors(adata, n_neighbors=n_neighbors)
                            st.success(f"✅ Υπολογίστηκαν γείτονες ({n_neighbors}) με βάση το PCA")

                            sc.tl.leiden(adata, resolution=resolution, flavor="igraph", directed=False, random_state=0)
                            st.success(f"✅ Εκτελέστηκε Clustering με Leiden (resolution={resolution})")
                            st.info("ℹ️ Το Leiden εντοπίζει ομάδες κυττάρων βάσει γειτνίασης.")

                            sc.tl.umap(adata, n_components=n_umap_components)
                            st.success("✅ UMAP ολοκληρώθηκε πάνω στο PCA")
                            st.info("ℹ️ Το UMAP δημιουργεί δισδιάστατη προβολή για εύκολη οπτικοποίηση.")

                    # Ενημέρωση session state
                    st.session_state["adata_analysis"] = adata.copy()

                    if "leiden" in adata.obs.columns:
                        st.success("✅ Τα clusters αποθηκεύτηκαν στο `adata.obs['leiden']`")
                        n_clusters = len(adata.obs["leiden"].unique())
                        st.info(f"📦 Βρέθηκαν {n_clusters} διαφορετικά clusters με resolution = {resolution}")


                    status.update(label="✅ Ανάλυση ολοκληρώθηκε", state="complete")

                except Exception as e:
                    st.error(f"❌ Σφάλμα: {str(e)}")
                    status.update(label="❌ Αποτυχία κατά την ανάλυση", state="error")

    
    # UMAP Προβολή
    if "adata_analysis" in st.session_state:
        st.info("ℹ️ Μπορείς να αλλάξεις χρωματισμό χωρίς να ξανατρέξεις την ανάλυση.")
        adata = st.session_state["adata_analysis"]
        color_options = ["batch", "celltype"] if "celltype" in adata.obs.columns else ["batch"]
        selected_color = st.selectbox("🎨 Χρωματισμός κατά", options=color_options, key="umap_color_by")

        # Προβολή αρχικού UMAP (μετά το preprocessing)
        if "adata_pre_umap" in st.session_state:
            st.subheader("🔍 Σύγκριση UMAP: Πριν & Μετά την Ανάλυση")
            col1, col2 = st.columns(2)

            with col1:
                st.markdown("#### 🕒 Πριν (μετά το Preprocessing)")
                adata_pre_umap = st.session_state["adata_pre_umap"]
                try:
                    if projection == "3D":
                        coords = adata_pre_umap.obsm["X_umap"]
                        labels = adata_pre_umap.obs[selected_color].astype(str).values

                        fig_pre = go.Figure()
                        for label in np.unique(labels):
                            mask = labels == label
                            fig_pre.add_trace(go.Scatter3d(
                                x=coords[mask, 0],
                                y=coords[mask, 1],
                                z=coords[mask, 2],
                                mode="markers",
                                name=str(label),
                                marker=dict(size=4),
                                hoverinfo="text"
                            ))
                        fig_pre.update_layout(
                            title="UMAP 3D (Before)",
                            scene=dict(xaxis_title="UMAP1", yaxis_title="UMAP2", zaxis_title="UMAP3"),
                            margin=dict(l=0, r=0, b=0, t=30),
                            height=600
                        )
                        st.plotly_chart(fig_pre, use_container_width=True)

                    else:
                        fig_pre = sc.pl.umap(adata_pre_umap, color=selected_color, return_fig=True, projection="2d")
                        st.pyplot(fig_pre)

                except Exception as e:
                    st.error(f"❌ Σφάλμα στην προβολή UMAP πριν την ανάλυση: {str(e)}")


        with col2:
            st.markdown("#### 🧭 Μετά (τελικό αποτέλεσμα)")

            try:
                if projection == "3D":
                    # Αν δεν έχει 3 διαστάσεις, επαναυπολογίζουμε UMAP
                    if adata.obsm["X_umap"].shape[1] < 3:
                        sc.tl.umap(adata, n_components=3)

                    coords = adata.obsm["X_umap"]
                    labels = adata.obs[selected_color].astype(str).values

                    fig = go.Figure()
                    for label in np.unique(labels):
                        mask = labels == label
                        fig.add_trace(go.Scatter3d(
                            x=coords[mask, 0],
                            y=coords[mask, 1],
                            z=coords[mask, 2],
                            mode="markers",
                            name=str(label),
                            marker=dict(size=4),
                            hoverinfo="text"
                        ))

                    fig.update_layout(
                        title=f"UMAP 3D ({selected_color})",
                        scene=dict(
                            xaxis_title="UMAP1",
                            yaxis_title="UMAP2",
                            zaxis_title="UMAP3"
                        ),
                        margin=dict(l=0, r=0, b=0, t=30),
                        height=600
                    )

                    st.plotly_chart(fig, use_container_width=True)

                else:
                    fig = sc.pl.umap(adata, color=selected_color, projection="2d", return_fig=True)
                    st.pyplot(fig)

            except Exception as e:
                st.error(f"❌ Σφάλμα στην προβολή {projection} UMAP: {str(e)}")


    st.markdown("---")
    st.info("➡️ Κάνε scroll επάνω και επίλεξε το επόμενο tab για να συνεχίσεις με την γονιδιακή ανάλυση.")

# Tab 4: Γονιδιακή Ανάλυση
with tab4:
    st.header("🧬 Γονιδιακή Ανάλυση")
    if "adata" not in st.session_state:
        st.warning("❗ Αρχικά ανέβασε αρχείο `.h5ad` στο πρώτο tab.")
    else:
        if "adata_analysis" not in st.session_state:
            st.warning("❗ Ολοκλήρωσε πρώτα την ανάλυση PCA/Clustering.")
        else:
            adata = st.session_state["adata_analysis"]

            # 1. Marker Genes
            st.subheader("🧪 Marker Genes")
            with st.expander("➕ Επιλογές Marker Genes", expanded=True):
                groupby_col = st.selectbox("🧩 Επιλέξτε πεδίο clustering", options=adata.obs.columns, key="mg_groupby")
                method = st.selectbox("⚙️ Μέθοδος κατάταξης", options=["wilcoxon", "logreg", "t-test"], key="mg_method")
                vis_type = st.radio("🎨 Επιλογή οπτικοποίησης", ["dotplot", "heatmap", "violin"], horizontal=True, key="vis_type_marker")
                run_marker = st.button("🔍 Εκτέλεση Marker Genes", key="run_marker_button")

            if run_marker:
                with st.status("🔄 Ανάλυση Marker Genes...", expanded=True) as status:
                    try:
                        if not pd.api.types.is_categorical_dtype(adata.obs[groupby_col]):
                            adata.obs[groupby_col] = adata.obs[groupby_col].astype("category")

                        group_sizes = adata.obs[groupby_col].value_counts()
                        if any(group_sizes < 2):
                            st.warning("⚠️ Κάποιες ομάδες έχουν λιγότερα από 2 δείγματα.")
                            status.update(label="❌ Ανάλυση ακυρώθηκε", state="error")
                            st.stop()

                        with st.spinner("🔍 Εκτελείται Marker Genes..."):
                            sc.tl.rank_genes_groups(adata, groupby=groupby_col, method=method, use_raw=False)
                            st.success("✅ Η ανάλυση ολοκληρώθηκε")

                        result = adata.uns["rank_genes_groups"]
                        groups = result["names"].dtype.names

                        st.session_state["marker_result"] = result
                        st.session_state["marker_groups"] = groups
                        st.session_state["marker_groupby_used"] = groupby_col
                        st.session_state["marker_raw"] = adata.copy()

                        status.update(label="✅ Marker Genes ολοκληρώθηκε", state="complete")

                    except Exception as e:
                        st.error(f"❌ Σφάλμα: {str(e)}")
                        status.update(label="❌ Αποτυχία Marker Genes", state="error")

            if (
                "marker_result" in st.session_state and
                "marker_groups" in st.session_state and
                st.session_state.get("marker_groupby_used") == st.session_state.get("mg_groupby")
            ):
                result = st.session_state["marker_result"]
                groups = st.session_state["marker_groups"]

                # 🔍 Οπτικοποιήσεις πρώτα
                groupby_col = st.session_state["mg_groupby"]
                adata = st.session_state["adata_analysis"]
                if not pd.api.types.is_categorical_dtype(adata.obs[groupby_col]):
                    adata.obs[groupby_col] = adata.obs[groupby_col].astype("category")

                try:
                    if st.session_state["vis_type_marker"] == "heatmap":
                        sc.tl.dendrogram(adata, groupby=groupby_col)
                        sc.pl.rank_genes_groups_heatmap(adata, groups=groups, n_genes=5, show=False)
                    elif st.session_state["vis_type_marker"] == "dotplot":
                        sc.pl.rank_genes_groups_dotplot(adata, groups=groups, n_genes=5, show=False)
                    elif st.session_state["vis_type_marker"] == "violin":
                        sc.pl.rank_genes_groups_violin(adata, groups=groups, n_genes=5, show=False)
                    st.pyplot(plt.gcf())
                except Exception as e:
                    st.warning("⚠️ Πάτησε ξανά 'Εκτέλεση Marker Genes' για σωστή εμφάνιση των αποτελεσμάτων.")
                    st.error(f"❌ Σφάλμα οπτικοποίησης: {str(e)}")

                # 👉 Data κάτω από το plot
                genes_per_page = st.selectbox("🔢 Αποτελέσματα ανά σελίδα", options=[10, 20, 50], index=0, key="marker_page_size")
                n_genes_total = len(result["names"][groups[0]])
                total_pages = (n_genes_total - 1) // genes_per_page + 1
                page = st.number_input("📄 Επιλογή σελίδας", min_value=1, max_value=total_pages, step=1, key="marker_page_num")

                start = (page - 1) * genes_per_page
                end = min(start + genes_per_page, n_genes_total)
                paginated_df = pd.DataFrame({group: result["names"][group][start:end] for group in groups})
                st.dataframe(paginated_df, use_container_width=True)

            elif "marker_result" in st.session_state:
                st.warning("⚠️ Το πεδίο clustering άλλαξε. Πάτησε ξανά 'Εκτέλεση Marker Genes' για σωστή εμφάνιση.")

            # 2. DEG Analysis
            st.subheader("📊 Ανάλυση Διαφορικής Έκφρασης (DEG)")
            with st.expander("➕ Επιλογές DEG", expanded=False):
                group_col = st.selectbox("🎯 Πεδίο σύγκρισης (group)", options=adata.obs.columns, key="deg_group")
                groups = adata.obs[group_col].unique().tolist()
                group = st.selectbox("🔬 Ομάδα ενδιαφέροντος", options=groups, key="deg_target")
                reference = st.selectbox("🆚 Ομάδα σύγκρισης (control)", options=[g for g in groups if g != group], key="deg_control")
                top_n = st.slider("🔝 Top UP/DOWN genes", min_value=5, max_value=100, value=20, step=5, key="deg_slider")
                run_deg = st.button("📊 Εκτέλεση DEG Ανάλυσης", key="run_deg_button")

            if run_deg:
                with st.status("🔄 Ανάλυση Διαφορικής Έκφρασης...", expanded=True) as status:
                    try:
                        if not pd.api.types.is_categorical_dtype(adata.obs[group_col]):
                            adata.obs[group_col] = adata.obs[group_col].astype("category")

                        group_sizes = adata.obs[group_col].value_counts()
                        if group_sizes[group] < 2 or group_sizes[reference] < 2:
                            st.warning("⚠️ Η ομάδα ενδιαφέροντος ή η ομάδα σύγκρισης έχει λιγότερα από 2 δείγματα.")
                            status.update(label="❌ Ανάλυση ακυρώθηκε", state="error")
                            st.stop()

                        with st.spinner("📊 Υπολογισμός DEG..."):
                            sc.tl.rank_genes_groups(
                                adata, groupby=group_col, groups=[group], reference=reference,
                                method="wilcoxon", use_raw=False
                            )
                            st.success("✅ DEG Analysis ολοκληρώθηκε")

                        deg = adata.uns["rank_genes_groups"]
                        if group in deg["names"].dtype.names:
                            degs_df = pd.DataFrame({
                                "genes": deg["names"][group],
                                "pvals": deg["pvals"][group],
                                "pvals_adj": deg["pvals_adj"][group],
                                "logfoldchanges": deg["logfoldchanges"][group],
                            })

                            degs_df["pvals"] = degs_df["pvals"].replace(0, np.nextafter(0, 1))
                            degs_df["neg_log10_pval"] = -np.log10(degs_df["pvals"])
                            degs_df["diffexpressed"] = "NS"
                            degs_df.loc[(degs_df["logfoldchanges"] > 1) & (degs_df["pvals"] < 0.05), "diffexpressed"] = "UP"
                            degs_df.loc[(degs_df["logfoldchanges"] < -1) & (degs_df["pvals"] < 0.05), "diffexpressed"] = "DOWN"

                            top_up = degs_df[degs_df["diffexpressed"] == "UP"].nlargest(top_n, "neg_log10_pval")
                            top_down = degs_df[degs_df["diffexpressed"] == "DOWN"].nlargest(top_n, "neg_log10_pval")

                            # 🟢 Πρώτα Volcano plot
                            st.subheader("🌋 Volcano Plot")
                            fig, ax = plt.subplots(figsize=(10, 6))
                            sns.scatterplot(
                                data=degs_df,
                                x="logfoldchanges",
                                y="neg_log10_pval",
                                hue="diffexpressed",
                                palette={"UP": "#bb0c00", "DOWN": "#00AFBB", "NS": "gray"},
                                alpha=0.7,
                                ax=ax,
                                edgecolor=None
                            )
                            ax.axhline(y=-np.log10(0.05), color='gray', linestyle='dashed')
                            ax.axvline(x=1, color='gray', linestyle='dashed')
                            ax.axvline(x=-1, color='gray', linestyle='dashed')
                            ax.set_xlabel("log2 Fold Change")
                            ax.set_ylabel("-log10 p-value")
                            ax.set_title(f"Volcano Plot: {group} vs {reference}")
                            st.pyplot(fig)

                            # 👉 Πίνακες μετά το plot
                            col1, col2 = st.columns(2)
                            col1.dataframe(top_up[["genes", "logfoldchanges", "pvals"]])
                            col2.dataframe(top_down[["genes", "logfoldchanges", "pvals"]])

                            st.session_state["deg_group_col_used"] = group_col
                            st.session_state["deg_group_used"] = group

                            status.update(label="✅ DEG Ανάλυση ολοκληρώθηκε", state="complete")

                            st.session_state["deg_group_selected"] = group
                            st.session_state["deg_top_n"] = top_n
                        else:
                            st.warning("⚠️ Δεν βρέθηκαν αποτελέσματα για την επιλεγμένη ομάδα.")
                    except Exception as e:
                        st.warning("⚠️ Πάτησε ξανά '📊 Εκτέλεση DEG Ανάλυσης' για σωστή εμφάνιση των αποτελεσμάτων.")
                        st.error(f"❌ Σφάλμα: {str(e)}")
                        status.update(label="❌ Αποτυχία DEG", state="error")

    st.markdown("---")
    st.info("➡️ Κάνε scroll επάνω και επίλεξε το επόμενο tab για να συνεχίσεις με την εξαγωγή αποτελεσμάτων.")


# Tab 5: Εξαγωγή Αποτελεσμάτων
with tab5:
    st.header("💾 Λήψεις Αποτελεσμάτων")

    # Λήψη preprocessed δεδομένων
    st.subheader("🧪 Preprocessed Δεδομένα")
    if "adata_pre" in st.session_state:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
            tmp_path = tmp.name
        st.session_state["adata_pre"].write(tmp_path)

        with open(tmp_path, "rb") as f:
            data = f.read()

        st.download_button(
            label="📥 Κατέβασε το Preprocessed αρχείο (.h5ad)",
            data=data,
            file_name="preprocessed_data.h5ad",
            mime="application/octet-stream"
        )

        os.remove(tmp_path)


    # DEG αποτελέσματα
    if "adata" not in st.session_state or "adata_analysis" not in st.session_state:
        st.warning("❗ Χρειάζεται να έχει προηγηθεί γονιδιακή ανάλυση.")
    elif "deg_group_selected" not in st.session_state or "deg_top_n" not in st.session_state:
        st.warning("❗ Χρειάζεται να έχει εκτελεστεί DEG στο Tab 4.")
    else:
        adata = st.session_state["adata_analysis"]
        group = st.session_state["deg_group_selected"]
        top_n = st.session_state["deg_top_n"]
        group_col = st.session_state["deg_group_col_used"]

        deg = adata.uns["rank_genes_groups"]

        if (
            st.session_state.get("deg_group_col_used") != group_col or
            st.session_state.get("deg_group_used") != group
        ):
            st.warning("⚠️ Το πεδίο σύγκρισης ή η ομάδα ενδιαφέροντος έχει αλλάξει.")
            st.stop()

        if group not in deg["names"].dtype.names:
            st.warning(f"⚠️ Δεν υπάρχουν αποτελέσματα για την ομάδα '{group}'.")
            st.stop()

        degs_df = pd.DataFrame({
            "genes": deg["names"][group],
            "pvals": deg["pvals"][group],
            "pvals_adj": deg["pvals_adj"][group],
            "logfoldchanges": deg["logfoldchanges"][group],
        })

        degs_df["pvals"] = degs_df["pvals"].replace(0, np.nextafter(0, 1))
        degs_df["neg_log10_pval"] = -np.log10(degs_df["pvals"])
        degs_df["diffexpressed"] = "NS"
        degs_df.loc[(degs_df["logfoldchanges"] > 1) & (degs_df["pvals"] < 0.05), "diffexpressed"] = "UP"
        degs_df.loc[(degs_df["logfoldchanges"] < -1) & (degs_df["pvals"] < 0.05), "diffexpressed"] = "DOWN"

        top_up = degs_df[degs_df["diffexpressed"] == "UP"].nlargest(top_n, "neg_log10_pval")
        top_down = degs_df[degs_df["diffexpressed"] == "DOWN"].nlargest(top_n, "neg_log10_pval")

        def to_excel(df):
            output = BytesIO()
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                df.to_excel(writer, index=False, sheet_name='DEGs')
            return output.getvalue()

        fig, ax = plt.subplots(figsize=(10, 6))
        sns.scatterplot(
            data=degs_df, x="logfoldchanges", y="neg_log10_pval",
            hue="diffexpressed",
            palette={"UP": "#bb0c00", "DOWN": "#00AFBB", "NS": "gray"},
            alpha=0.7, ax=ax
        )
        ax.axhline(y=-np.log10(0.05), color='gray', linestyle='dashed')
        ax.axvline(x=1, color='gray', linestyle='dashed')
        ax.axvline(x=-1, color='gray', linestyle='dashed')
        ax.set_title(f"Volcano Plot: {group}")
        buf = BytesIO()
        fig.savefig(buf, format="png", dpi=300)

        st.subheader(f"📂 Λήψεις για ομάδα: `{group}`  (Top {top_n} genes)")
        col1, col2 = st.columns(2)
        with col1:
            st.download_button("📥 Όλα τα DEGs (.csv)", data=degs_df.to_csv(index=False), mime="text/csv", file_name=f"DEGs_{group}.csv")
            st.download_button("📥 Top UP (.csv)", data=top_up.to_csv(index=False), mime="text/csv", file_name=f"Top{top_n}_UP_{group}.csv")
            st.download_button("📥 Top DOWN (.csv)", data=top_down.to_csv(index=False), mime="text/csv", file_name=f"Top{top_n}_DOWN_{group}.csv")
            st.download_button("📷 Volcano (.png)", data=buf.getvalue(), file_name=f"Volcano_{group}.png", mime="image/png")
        with col2:
            st.download_button("📥 Όλα τα DEGs (.xlsx)", data=to_excel(degs_df), mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", file_name=f"DEGs_{group}.xlsx")
            st.download_button("📥 Top UP (.xlsx)", data=to_excel(top_up), mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", file_name=f"Top{top_n}_UP_{group}.xlsx")
            st.download_button("📥 Top DOWN (.xlsx)", data=to_excel(top_down), mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", file_name=f"Top{top_n}_DOWN_{group}.xlsx")

        # --- UMAP εικόνες ---
        st.subheader("🧭 Λήψεις UMAP Διαγραμμάτων")
        try:
            fig_pre = sc.pl.umap(st.session_state["adata_pre_umap"], color=st.session_state["umap_color_by"], return_fig=True)
            buf_pre = BytesIO()
            fig_pre.savefig(buf_pre, format="png", dpi=300)
            st.download_button("📷 UMAP Πριν (.png)", data=buf_pre.getvalue(), file_name="umap_before.png", mime="image/png")
        except:
            st.warning("⚠️ Δεν υπάρχει διαθέσιμο UMAP Πριν.")

        try:
            fig_after = sc.pl.umap(adata, color=st.session_state["umap_color_by"], return_fig=True)
            buf_after = BytesIO()
            fig_after.savefig(buf_after, format="png", dpi=300)
            st.download_button("📷 UMAP Μετά (.png)", data=buf_after.getvalue(), file_name="umap_after.png", mime="image/png")
        except:
            st.warning("⚠️ Δεν υπάρχει διαθέσιμο UMAP Μετά.")

        # --- Οπτικοποιήσεις Marker Genes ---
        st.subheader("🎨 Οπτικοποιήσεις Marker Genes")
        if (
            "marker_groups" in st.session_state and
            "marker_raw" in st.session_state and
            "marker_groupby_used" in st.session_state
        ):
            try:
                adata = st.session_state["marker_raw"].copy()
                marker_groups = st.session_state["marker_groups"]
                groupby_col = st.session_state["marker_groupby_used"]

                if groupby_col not in adata.obs.columns:
                    st.warning(f"⚠️ Το πεδίο '{groupby_col}' δεν υπάρχει στο AnnData για Marker Genes.")
                    st.stop()

                if not pd.api.types.is_categorical_dtype(adata.obs[groupby_col]):
                    adata.obs[groupby_col] = adata.obs[groupby_col].astype("category")

                # Heatmap (μόνο αν υπάρχουν >2 groups)
                if len(marker_groups) > 2:
                    sc.tl.dendrogram(adata, groupby=groupby_col)
                    sc.pl.rank_genes_groups_heatmap(adata, groups=marker_groups, n_genes=5, show=False)
                    buf = BytesIO()
                    plt.savefig(buf, format="png", dpi=300)
                    st.download_button("📷 Heatmap (.png)", data=buf.getvalue(), file_name="marker_heatmap.png", mime="image/png")
                    plt.clf()

                # Dotplot
                sc.pl.rank_genes_groups_dotplot(adata, groups=marker_groups, n_genes=5, show=False)
                buf = BytesIO()
                plt.savefig(buf, format="png", dpi=300)
                st.download_button("📷 Dotplot (.png)", data=buf.getvalue(), file_name="marker_dotplot.png", mime="image/png")
                plt.clf()

                # Violin
                sc.pl.rank_genes_groups_violin(adata, groups=marker_groups, n_genes=5, show=False)
                buf = BytesIO()
                plt.savefig(buf, format="png", dpi=300)
                st.download_button("📷 Violin (.png)", data=buf.getvalue(), file_name="marker_violin.png", mime="image/png")
                plt.clf()

            except Exception as e:
                st.warning(f"⚠️ Σφάλμα στην εξαγωγή Marker Genes plots: {str(e)}")



    st.markdown("---")
    st.info("➡️ Κάνε scroll επάνω και επίλεξε το επόμενο tab για να συνεχίσεις με τις πληροφορίες για την ομάδα.")


with tab6:
    st.header("👥 Ομάδα Εργασίας")

    col_logo, col_info = st.columns([1, 3])
    with col_logo:
        st.image("images/ionian_university.png", width=120)
    with col_info:
        st.markdown("#### 🎓 Πληροφορίες Εφαρμογής")
        st.markdown("""
        Η εφαρμογή αναπτύχθηκε στο πλαίσιο της ακαδημαϊκής εργασίας για το μάθημα 
        **Τεχνολογίες Λογισμικού** στο **6ο εξάμηνο** σπουδών, 
        στο **Τμήμα Πληροφορικής του Ιονίου Πανεπιστημίου**, το 2025.
        """)
    st.divider()

    st.markdown("#### 🤝 Συνεργασία Ομάδας")
    st.markdown("""
    Όλα τα μέλη συμμετείχαν ενεργά στην υλοποίηση της εφαρμογής, τόσο στον προγραμματισμό όσο και στη συγγραφή της τεκμηρίωσης. 
    Παρακάτω αναφέρεται η βασική συνεισφορά καθενός ανά θεματική περιοχή:
    """)

    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("##### 👩‍💻 Ευτυχία Φίλιου")
        st.markdown("- Α.Μ.: `Π2019202`")
        st.markdown("- 📁 Tab 1 (Δεδομένα) & ⚙️ Tab 2 (Προεπεξεργασία)")

    with col2:
        st.markdown("##### 👨‍💻 Χρήστος Σπυρίδων Καρύδης")
        st.markdown("- Α.Μ.: `INF2022076`")
        st.markdown("- 📊 Tab 3 (Ανάλυση) & 🧬 Tab 4 (Γονιδιακή Ανάλυση)")

    with col3:
        st.markdown("##### 👩‍💻 Ευαγγελία Φώτη")
        st.markdown("- Α.Μ.: `INF2022224`")
        st.markdown("- 📈 Tab 5 (Αποτελέσματα) & 💾 Tab 6 (Ομάδα)")

    st.divider()

    st.markdown("#### 🖼️ Ομαδική Φωτογραφία")
    st.image("images/team.jpg", caption="Η ομάδα ανάπτυξης της εφαρμογής")

    st.divider()

    st.markdown("#### 🔗 Αποθετήριο Κώδικα")
    st.markdown("Το πλήρες αποθετήριο του έργου είναι διαθέσιμο στο GitHub:")
    st.markdown("[📂 Μετάβαση στο GitHub](https://github.com/chriskarydis/scRNA_seq_Pipeline)")

    st.divider()

