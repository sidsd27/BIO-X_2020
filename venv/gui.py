from PyQt5 import uic
from PyQt5.QtWidgets import QApplication



Form, Window = uic.loadUiType("/Users/sidsd27/Desktop/BIO-X_gui.ui")

app = QApplication([])
window = Window()
form = Form()
form.setupUi(window)
window.show()
app.exec_()