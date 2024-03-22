from pymol import cmd


def refresh_dropdown(dropdown_to_refresh):
    list_pymol_objects = cmd.get_names('all')
    dropdown_to_refresh.clear()
    dropdown_to_refresh.addItems(list_pymol_objects)
    if len(list_pymol_objects) > 0:
        dropdown_to_refresh.setCurrentText(list_pymol_objects[0])