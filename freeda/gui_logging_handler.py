"""

adapted from https://beenje.github.io/blog/posts/logging-to-a-tkinter-scrolledtext-widget/

Handles logging into Tkinter ScrollText widget

"""

import tkinter
import logging


class TextHandler(logging.Handler):
    """Logs messages into the scroll widget of the GUI"""

    def __init__(self, text):
        logging.Handler.__init__(self)
        self.text = text

    def emit(self, record):
        msg = self.format(record)

        def append():
            self.text.configure(state='normal')
            self.text.insert(tkinter.END, msg + '\n')
            self.text.configure(state='disabled')
            self.text.yview(tkinter.END)
        self.text.after(0, append)

