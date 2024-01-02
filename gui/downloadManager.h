#ifndef DOWNLOADMANAGER_H
#define DOWNLOADMANAGER_H

#include <QObject>
#include <QByteArray>
#include <QNetworkAccessManager>
#include <QNetworkRequest>
#include <QNetworkReply>

class downloadManager : public QObject
{
    Q_OBJECT
public:

    explicit downloadManager(QUrl imageUrl, QObject *parent = 0);
    virtual ~downloadManager();
    QByteArray downloadedData() const;
    bool saveToDisk(const QString &filename, QIODevice *data);
    QString saveFileName(const QUrl &url);

signals:

    void downloaded();

private slots:

    void fileDownloaded(QNetworkReply* pReply);

private:

    QNetworkAccessManager m_WebCtrl;
    QByteArray m_DownloadedData;

};


#endif // DOWNLOADMANAGER_H
