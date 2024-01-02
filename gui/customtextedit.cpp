// Copyright (C) 2016 The Qt Company Ltd.
// SPDX-License-Identifier: LicenseRef-Qt-Commercial OR BSD-3-Clause

#include "gui/customtextedit.h"
#include <QCompleter>
#include <QKeyEvent>
#include <QAbstractItemView>
#include <QtDebug>
#include <QApplication>
#include <QModelIndex>
#include <QAbstractItemModel>
#include <QScrollBar>

//! [0]
CustomTextEdit::CustomTextEdit(QWidget *parent): QPlainTextEdit(parent)
{
    //setPlainText(tr("This CustomTextEdit provides autocompletions for words that have more than" " 3 characters. You can trigger autocompletion using ") +
    //              QKeySequence("Ctrl+E").toString(QKeySequence::NativeText));
}
//! [0]

//! [1]
CustomTextEdit::~CustomTextEdit()
{
}
//! [1]

//! [2]
void CustomTextEdit::setCompleter(QCompleter *completer)
{
    if (c)
        c->disconnect(this);

    c = completer;

    if (!c)
        return;

    c->setWidget(this);
    c->setCompletionMode(QCompleter::PopupCompletion);
    c->setCaseSensitivity(Qt::CaseSensitive);
    QObject::connect(c, QOverload<const QString &>::of(&QCompleter::activated), this, &CustomTextEdit::insertCompletion);
}
//! [2]

//! [3]
QCompleter *CustomTextEdit::completer() const
{
    return c;
}
//! [3]

//! [4]
void CustomTextEdit::insertCompletion(const QString &completion)
{
    if (c->widget() != this)
        return;
    QTextCursor tc = textCursor();
    int extra = completion.length() - c->completionPrefix().length();
    tc.movePosition(QTextCursor::Left);
    tc.movePosition(QTextCursor::EndOfWord);
    tc.insertText(completion.right(extra));
    setTextCursor(tc);

    //QTextStream(stdout) << "completion:" << completion << " extra:" << extra << " completion.length():" << completion.length() << " c->completionPrefix().length():" << c->completionPrefix().length() << "\n";

}
//! [4]

//! [5]
QString CustomTextEdit::textUnderCursor() const
{
    QTextCursor tc = textCursor();
    tc.select(QTextCursor::WordUnderCursor);
    //QTextStream(stdout) << "QTextCursor::WordUnderCursor:" << QTextCursor::WordUnderCursor << " tc.selectedText():" << tc.selectedText()<< "\n";

    return tc.selectedText();
}
//! [5]

//! [6]
void CustomTextEdit::focusInEvent(QFocusEvent *e)
{
    if (c)
        c->setWidget(this);
    QPlainTextEdit::focusInEvent(e);
}
//! [6]

//! [7]
void CustomTextEdit::keyPressEvent(QKeyEvent *e)
{
    if (c && c->popup()->isVisible()) {
        // The following keys are forwarded by the completer to the widget
        switch (e->key()) {
        case Qt::Key_Enter:
        case Qt::Key_Return:
        case Qt::Key_Escape:
        case Qt::Key_Tab:
        case Qt::Key_Backtab:
            e->ignore();
            return; // let the completer do default behavior
        default:
            break;
        }
    }

    const bool isShortcut = (e->modifiers().testFlag(Qt::ControlModifier) && e->key() == Qt::Key_E); // CTRL+E
    if (!c || !isShortcut) // do not process the shortcut when we have a completer
        QPlainTextEdit::keyPressEvent(e);
    //! [7]

    //! [8]
    const bool ctrlOrShift = e->modifiers().testFlag(Qt::ControlModifier) ||
            e->modifiers().testFlag(Qt::ShiftModifier);
    if (!c || (ctrlOrShift && e->text().isEmpty()))
        return;

    static QString eow("~!@#$%^&*()+{}|:\"<>?,./;'[]\\-="); // end of word, i removed "_"
    const bool hasModifier = (e->modifiers() != Qt::NoModifier) && !ctrlOrShift;
    QString completionPrefix = textUnderCursor();

    if (!isShortcut && (hasModifier || e->text().isEmpty()|| completionPrefix.length() < 1
                        || eow.contains(e->text().right(1)))) {
        c->popup()->hide();
        return;
    }

    if (completionPrefix != c->completionPrefix()) {
        c->setCompletionPrefix(completionPrefix);
        c->popup()->setCurrentIndex(c->completionModel()->index(0, 0));
    }
    QRect cr = cursorRect();
    cr.setWidth(c->popup()->sizeHintForColumn(0) + c->popup()->verticalScrollBar()->sizeHint().width());
    c->complete(cr); // popup it up!
}
//! [8]

