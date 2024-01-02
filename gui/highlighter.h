#ifndef HIGHLIGHTER_H
#define HIGHLIGHTER_H

#include "qsyntaxhighlighter.h"
#include "QRegularExpression"

class Highlighter : public QSyntaxHighlighter
{
    Q_OBJECT

public:
    Highlighter(QTextDocument *parent = 0);

protected:
    void highlightBlock(const QString &text) override;

private:
    struct HighlightingRule
    {
        QRegularExpression pattern;
        QTextCharFormat format;
    };
    QVector<HighlightingRule> highlightingRules;

    QRegularExpression commentStartExpression;
    QRegularExpression commentEndExpression;

    QTextCharFormat keywordFormat;
    QTextCharFormat classFormat;
    QTextCharFormat singleLineCommentFormat;
    QTextCharFormat multiLineCommentFormat;
    QTextCharFormat quotationFormat;
    QTextCharFormat functionFormat;

    QTextCharFormat DirectoryCommandFormat;
    QTextCharFormat CommandFormat;
    QTextCharFormat NumPointsFormat;
    QTextCharFormat UnitsFormat;
    QTextCharFormat ConfigFormat;
    QTextCharFormat ConstValuesFormat;

};

#endif // HIGHLIGHTER_H
